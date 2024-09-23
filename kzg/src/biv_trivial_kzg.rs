use ark_ec::{
    pairing::{Pairing, PairingOutput},
    scalar_mul::variable_base::VariableBaseMSM,
    CurveGroup, Group,
};
use ark_ff::{Field, One, UniformRand, Zero};
use ark_poly::polynomial::{
    univariate::DensePolynomial as UnivariatePolynomial, DenseUVPolynomial, Polynomial,
};

use crate::trivial_kzg::structured_generators_scalar_power;

use ark_std::{end_timer, start_timer};
use std::{marker::PhantomData, ops::Mul};

use ark_std::rand::Rng;
use digest::Digest;

use crate::Error;

pub struct BivariatePolynomial<F: Field> {
    x_polynomials: Vec<UnivariatePolynomial<F>>,
}

// We want \sum_i f_i(X) Y^{i-1}
impl<F: Field> BivariatePolynomial<F> {
    pub fn evaluate(&self, point: &(F, F)) -> F {
        let (x, y) = point;
        let mut point_y_powers = vec![];
        let mut cur = F::one();
        for _ in 0..(self.x_polynomials.len()) {
            point_y_powers.push(cur);
            cur *= y;
        }
        point_y_powers
            .iter()
            .zip(&self.x_polynomials)
            .map(|(y_power, x_polynomial)| y_power.clone() * x_polynomial.evaluate(&x))
            .sum()
    }
}

pub struct BivariatePolynomialCommitment<P: Pairing> {
    _pairing: PhantomData<P>,
}

impl<P: Pairing> BivariatePolynomialCommitment<P> {
    pub fn setup<R: Rng>(
        rng: &mut R,
        x_degree: usize,
        y_degree: usize,
    ) -> Result<(Vec<Vec<P::G1Affine>>), Error> {
        let alpha = <P::ScalarField>::rand(rng);
        let beta = <P::ScalarField>::rand(rng);
        let gamma = <P::ScalarField>::rand(rng);
        let g = <P::G1>::generator();
        let h = g.mul(beta);

        let k = g + h;
        let x_srs = <P as Pairing>::G1::normalize_batch(&structured_generators_scalar_power(
            x_degree + 1,
            &g,
            &alpha,
        ));
        let y_srs = <P as Pairing>::G1::normalize_batch(&structured_generators_scalar_power(
            x_degree + 1,
            &h,
            &gamma,
        ));


        let srs = SRS {
            g_alpha_powers: vec![g.clone()],
            h_beta_powers: structured_generators_scalar_power(2 * x_degree + 1, &h, &beta),
            g_beta: g * beta,
            h_alpha: h * alpha,
        };
        Ok((srs, kzg_srs))
    }

    pub fn commit(
        srs: &(SRS<P>, Vec<P::G1Affine>),
        bivariate_polynomial: &BivariatePolynomial<P::ScalarField>,
    ) -> Result<(PairingOutput<P>, Vec<P::G1>), Error> {
        let (ip_srs, kzg_srs) = srs;
        let (ck, _) = ip_srs.get_commitment_keys();
        assert!(ck.len() >= bivariate_polynomial.x_polynomials.len());

        // Create KZG commitments to Y polynomials
        let x_polynomial_coms = bivariate_polynomial
            .x_polynomials
            .iter()
            .chain(vec![UnivariatePolynomial::zero()].iter().cycle())
            .take(ck.len())
            .map(|x_polynomial| KZG::<P>::commit(kzg_srs, x_polynomial))
            .collect::<Result<Vec<P::G1>, Error>>()?;

        // Create AFGHO commitment to Y polynomial commitments
        Ok((
            AFGHOCommitmentG1::<P>::commit(&ck, &x_polynomial_coms)?,
            x_polynomial_coms,
        ))
    }

    pub fn open(
        srs: &(SRS<P>, Vec<P::G1Affine>),
        bivariate_polynomial: &BivariatePolynomial<P::ScalarField>,
        x_polynomial_comms: &Vec<P::G1>,
        point: &(P::ScalarField, P::ScalarField),
    ) -> Result<OpeningProof<P, D>, Error> {
        let (x, y) = point;
        let (ip_srs, kzg_srs) = srs;
        let (ck_1, _) = ip_srs.get_commitment_keys();
        assert!(ck_1.len() >= bivariate_polynomial.x_polynomials.len());

        let precomp_time = start_timer!(|| "Computing coefficients and KZG commitment");
        let mut powers_of_x = vec![];
        let mut cur = P::ScalarField::one();
        for _ in 0..(ck_1.len()) {
            powers_of_x.push(cur);
            cur *= x;
        }

        let coeffs = bivariate_polynomial
            .x_polynomials
            .iter()
            .chain(vec![UnivariatePolynomial::zero()].iter().cycle())
            .take(ck_1.len())
            .map(|x_polynomial| {
                let mut c = x_polynomial.coeffs.to_vec();
                c.resize(kzg_srs.len(), <P::ScalarField>::zero());
                c
            })
            .collect::<Vec<Vec<P::ScalarField>>>();
        let y_eval_coeffs = (0..kzg_srs.len())
            .map(|j| {
                (0..ck_1.len())
                    .map(|i| powers_of_x[i].clone() * &coeffs[i][j])
                    .sum()
            })
            .collect::<Vec<P::ScalarField>>();
        // Can unwrap because y_eval_coeffs.len() is guarnateed to be equal to kzg_srs.len()
        let y_eval_comm = P::G1::msm(kzg_srs, &y_eval_coeffs).unwrap();
        end_timer!(precomp_time);

        let ipa_time = start_timer!(|| "Computing IPA proof");
        let ip_proof =
            PolynomialEvaluationSecondTierIPA::<P, D>::prove_with_structured_scalar_message(
                &ip_srs,
                (x_polynomial_comms, &powers_of_x),
                (&ck_1, &HomomorphicPlaceholderValue),
            )?;
        end_timer!(ipa_time);
        let kzg_time = start_timer!(|| "Computing KZG opening proof");
        let kzg_proof = KZG::<P>::open(
            kzg_srs,
            &UnivariatePolynomial::from_coefficients_slice(&y_eval_coeffs),
            y,
        )?;
        end_timer!(kzg_time);

        Ok(OpeningProof {
            ip_proof,
            y_eval_comm,
            kzg_proof,
        })
    }

    pub fn verify(
        v_srs: &VerifierSRS<P>,
        com: &PairingOutput<P>,
        point: &(P::ScalarField, P::ScalarField),
        eval: &P::ScalarField,
        proof: &OpeningProof<P, D>,
    ) -> Result<bool, Error> {
        let (x, y) = point;
        let ip_proof_valid =
            PolynomialEvaluationSecondTierIPA::<P, D>::verify_with_structured_scalar_message(
                v_srs,
                &HomomorphicPlaceholderValue,
                (com, &IdentityOutput(vec![proof.y_eval_comm.clone()])),
                x,
                &proof.ip_proof,
            )?;
        let kzg_proof_valid =
            KZG::<P>::verify(v_srs, &proof.y_eval_comm, y, eval, &proof.kzg_proof)?;
        Ok(ip_proof_valid && kzg_proof_valid)
    }
}

pub struct UnivariatePolynomialCommitment<P: Pairing, D: Digest> {
    _pairing: PhantomData<P>,
    _digest: PhantomData<D>,
}

impl<P: Pairing, D: Digest> UnivariatePolynomialCommitment<P, D> {
    fn bivariate_degrees(univariate_degree: usize) -> (usize, usize) {
        //(((univariate_degree + 1) as f64).sqrt().ceil() as usize).next_power_of_two() - 1;
        let sqrt = (((univariate_degree + 1) as f64).sqrt().ceil() as usize).next_power_of_two();
        // Skew split between bivariate degrees to account for KZG being less expensive than MIPP
        let skew_factor = if sqrt >= 32 { 16_usize } else { sqrt / 2 };
        (sqrt / skew_factor - 1, sqrt * skew_factor - 1)
    }

    fn parse_bivariate_degrees_from_srs(srs: &(SRS<P>, Vec<P::G1Affine>)) -> (usize, usize) {
        let x_degree = (srs.0.h_beta_powers.len() - 1) / 2;
        let y_degree = srs.1.len() - 1;
        (x_degree, y_degree)
    }

    fn bivariate_form(
        bivariate_degrees: (usize, usize),
        polynomial: &UnivariatePolynomial<P::ScalarField>,
    ) -> BivariatePolynomial<P::ScalarField> {
        let (x_degree, y_degree) = bivariate_degrees;
        let default_zero = vec![P::ScalarField::zero()];
        let mut coeff_iter = polynomial
            .coeffs
            .iter()
            .chain(default_zero.iter().cycle())
            .take((x_degree + 1) * (y_degree + 1));

        let mut x_polynomials = Vec::new();
        for _ in 0..x_degree + 1 {
            let mut x_polynomial_coeffs = vec![];
            for _ in 0..y_degree + 1 {
                x_polynomial_coeffs.push(Clone::clone(coeff_iter.next().unwrap()))
            }
            x_polynomials.push(UnivariatePolynomial::from_coefficients_slice(
                &x_polynomial_coeffs,
            ));
        }
        BivariatePolynomial { x_polynomials }
    }

    pub fn setup<R: Rng>(rng: &mut R, degree: usize) -> Result<(SRS<P>, Vec<P::G1Affine>), Error> {
        let (x_degree, y_degree) = Self::bivariate_degrees(degree);
        BivariatePolynomialCommitment::<P, D>::setup(rng, x_degree, y_degree)
    }

    pub fn commit(
        srs: &(SRS<P>, Vec<P::G1Affine>),
        polynomial: &UnivariatePolynomial<P::ScalarField>,
    ) -> Result<(PairingOutput<P>, Vec<P::G1>), Error> {
        let bivariate_degrees = Self::parse_bivariate_degrees_from_srs(srs);
        BivariatePolynomialCommitment::<P, D>::commit(
            srs,
            &Self::bivariate_form(bivariate_degrees, polynomial),
        )
    }

    pub fn open(
        srs: &(SRS<P>, Vec<P::G1Affine>),
        polynomial: &UnivariatePolynomial<P::ScalarField>,
        x_polynomial_comms: &Vec<P::G1>,
        point: &P::ScalarField,
    ) -> Result<OpeningProof<P, D>, Error> {
        let (x_degree, y_degree) = Self::parse_bivariate_degrees_from_srs(srs);
        let y = point.clone();
        let x = point.pow(&vec![(y_degree + 1) as u64]);
        BivariatePolynomialCommitment::open(
            srs,
            &Self::bivariate_form((x_degree, y_degree), polynomial),
            x_polynomial_comms,
            &(x, y),
        )
    }

    pub fn verify(
        v_srs: &VerifierSRS<P>,
        max_degree: usize,
        com: &PairingOutput<P>,
        point: &P::ScalarField,
        eval: &P::ScalarField,
        proof: &OpeningProof<P, D>,
    ) -> Result<bool, Error> {
        let (_, y_degree) = Self::bivariate_degrees(max_degree);
        let y = point.clone();
        let x = y.pow(&vec![(y_degree + 1) as u64]);
        BivariatePolynomialCommitment::verify(v_srs, com, &(x, y), eval, proof)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_std::rand::{rngs::StdRng, SeedableRng};
    use blake2::Blake2b;

    const BIVARIATE_X_DEGREE: usize = 7;
    const BIVARIATE_Y_DEGREE: usize = 7;
    //const UNIVARIATE_DEGREE: usize = 56;
    const UNIVARIATE_DEGREE: usize = 65535;
    //const UNIVARIATE_DEGREE: usize = 1048575;

    type TestBivariatePolyCommitment = BivariatePolynomialCommitment<Bls12_381, Blake2b>;
    type TestUnivariatePolyCommitment = UnivariatePolynomialCommitment<Bls12_381, Blake2b>;

    #[test]
    fn bivariate_poly_commit_test() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let srs =
            TestBivariatePolyCommitment::setup(&mut rng, BIVARIATE_X_DEGREE, BIVARIATE_Y_DEGREE)
                .unwrap();
        let v_srs = srs.0.get_verifier_key();

        let mut x_polynomials = Vec::new();
        for _ in 0..BIVARIATE_X_DEGREE + 1 {
            let mut x_polynomial_coeffs = vec![];
            for _ in 0..BIVARIATE_Y_DEGREE + 1 {
                x_polynomial_coeffs.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
            }
            x_polynomials.push(UnivariatePolynomial::from_coefficients_slice(
                &x_polynomial_coeffs,
            ));
        }
        let bivariate_polynomial = BivariatePolynomial { x_polynomials };

        // Commit to polynomial
        let (com, x_polynomial_comms) =
            TestBivariatePolyCommitment::commit(&srs, &bivariate_polynomial).unwrap();

        // Evaluate at challenge point
        let point = (UniformRand::rand(&mut rng), UniformRand::rand(&mut rng));
        let eval_proof = TestBivariatePolyCommitment::open(
            &srs,
            &bivariate_polynomial,
            &x_polynomial_comms,
            &point,
        )
        .unwrap();
        let eval = bivariate_polynomial.evaluate(&point);

        // Verify proof
        assert!(
            TestBivariatePolyCommitment::verify(&v_srs, &com, &point, &eval, &eval_proof).unwrap()
        );
    }

    // `cargo test univariate_poly_commit_test --release --features print-trace -- --ignored --nocapture`
    #[ignore]
    #[test]
    fn univariate_poly_commit_test() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let srs = TestUnivariatePolyCommitment::setup(&mut rng, UNIVARIATE_DEGREE).unwrap();
        let v_srs = srs.0.get_verifier_key();

        let mut polynomial_coeffs = vec![];
        for _ in 0..UNIVARIATE_DEGREE + 1 {
            polynomial_coeffs.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
        }
        let polynomial = UnivariatePolynomial::from_coefficients_slice(&polynomial_coeffs);

        // Commit to polynomial
        let (com, x_polynomial_comms) =
            TestUnivariatePolyCommitment::commit(&srs, &polynomial).unwrap();

        // Evaluate at challenge point
        let point = UniformRand::rand(&mut rng);
        let eval_proof =
            TestUnivariatePolyCommitment::open(&srs, &polynomial, &x_polynomial_comms, &point)
                .unwrap();
        let eval = polynomial.evaluate(&point);

        // Verify proof
        assert!(TestUnivariatePolyCommitment::verify(
            &v_srs,
            UNIVARIATE_DEGREE,
            &com,
            &point,
            &eval,
            &eval_proof
        )
        .unwrap());
    }
}
