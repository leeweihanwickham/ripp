use ark_ec::{
    pairing::Pairing,
    CurveGroup, 
    Group,
    scalar_mul::variable_base::VariableBaseMSM,
    // scalar_mul::fixed_base::FixedBase,
};
use ark_ff::{One, Field, UniformRand, Zero};
use ark_poly::{polynomial::{
    univariate::DensePolynomial as UnivariatePolynomial, Polynomial,
}, DenseUVPolynomial};

use crate::trivial_kzg::structured_generators_scalar_power;

use std::marker::PhantomData;

use ark_std::rand::Rng;
// use digest::Digest;

use crate::Error;

pub struct VerifierSRS<P: Pairing> {
    pub g: P::G1,
    pub h: P::G2,
    pub h_alpha: P::G2,
    pub h_beta: P::G2
}

pub struct BivariatePolynomial<F: Field> {
    pub x_polynomials: Vec<UnivariatePolynomial<F>>,
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

pub struct BivariateBatchKZG<P: Pairing> {
    _pairing: PhantomData<P>,
}

impl<P: Pairing> BivariateBatchKZG<P> {
    pub fn setup<R: Rng>(
        rng: &mut R,
        x_degree: usize,
        y_degree: usize,
    ) -> Result<(Vec<Vec<P::G1Affine>>, VerifierSRS<P>), Error> {
        let alpha = <P::ScalarField>::rand(rng);
        let beta = <P::ScalarField>::rand(rng);
        let g = <P::G1>::generator();
        let h = <P::G2>::generator();

        let mut final_srs: Vec<Vec<P::G1Affine>> = Vec::new();
        let mut temp = g;
        for _ in 0..(y_degree + 1) {
            let temp_srs = structured_generators_scalar_power(
                x_degree + 1,
                &temp,
                &alpha,
            );
            final_srs.push( <P as Pairing>::G1::normalize_batch(&temp_srs));
            temp *= beta;
        }

        Ok((final_srs,
            VerifierSRS {
                g: g.clone(),
                h: h.clone(),
                h_alpha: h * alpha,
                h_beta: h * beta
        }))
    }

    pub fn commit(
        powers: &Vec<Vec<P::G1Affine>>,
        bivariate_polynomial: &BivariatePolynomial<P::ScalarField>,
    ) -> Result<P::G1, Error> {
        assert!(powers.len() == bivariate_polynomial.x_polynomials.len());
        assert!(powers[0].len() >= bivariate_polynomial.x_polynomials[0].degree() + 1);

        let mut extended_coeff: Vec<<P as Pairing>::ScalarField> = Vec::new();
        let mut extended_powers: Vec<<P as Pairing>::G1Affine> = Vec::new();

        for i in 0..powers.len() {
            let mut coeffs = bivariate_polynomial.x_polynomials[i].coeffs.to_vec();
            coeffs.resize(powers[0].len(), <P::ScalarField>::zero());
            extended_coeff.extend(&coeffs);

            extended_powers.extend(&powers[i]);
        }
        
        Ok(P::G1::msm(&extended_powers, &extended_coeff).unwrap())
    }

    pub fn open(
        powers: &Vec<Vec<P::G1Affine>>,
        bivariate_polynomial: &BivariatePolynomial<P::ScalarField>,
        point: &(P::ScalarField, P::ScalarField),
    ) -> Result<(P::G1, P::G1), Error> {
        // let (x_srs, y_srs) = (&powers.0[0], &powers.1);
        // assert_eq!(x_srs.len(), bivariate_polynomial.x_polynomials[0].len());
        // assert_eq!(y_srs.len(), bivariate_polynomial.x_polynomials.len());

        let y_srs: Vec<<P as Pairing>::G1Affine> = powers.iter()
            .filter_map(|row| row.get(0))
            .cloned()
            .collect();

        // assert_eq!(y_srs, &y_srs_test);

        // let f(x,y) = \sum_i f_i(x) y^{i-1}
        // let p1(x, y) = f(x, y) - f(z1, y) / (x - z1)
        // p2(x, y) = f(z1, y) - f(z1, z2) / (y - z2)

        // compute f_i(z1)
        let evals: Vec<P::ScalarField> = bivariate_polynomial.x_polynomials
            .iter()
            .map(|poly| poly.evaluate(&point.0))
            .collect();

        let mut extended_coeff: Vec<<P as Pairing>::ScalarField> = Vec::new();
        let mut extended_powers: Vec<<P as Pairing>::G1Affine> = Vec::new();

        for i in 0..y_srs.len() {
            let x_polynomial = bivariate_polynomial.x_polynomials[i].clone() +
                UnivariatePolynomial::from_coefficients_vec(vec![-evals[i]]);
            let quotient_polynomial_x = &x_polynomial
                / &UnivariatePolynomial::from_coefficients_vec(vec![
                    -point.0.clone(),
                    P::ScalarField::one()
                ]);
            let mut quotient_coeffs_x = quotient_polynomial_x.coeffs.to_vec();
            quotient_coeffs_x.resize(powers[0].len(), <P::ScalarField>::zero());
            
            extended_coeff.extend(&quotient_coeffs_x);
            extended_powers.extend(&powers[i]);
        }
        // let quotient_polynomial_y = BivariatePolynomial { x_polynomials };

        let polynomial_z1_y = UnivariatePolynomial::from_coefficients_vec(evals);
        let quotient_polynomial_y = &polynomial_z1_y 
            / &UnivariatePolynomial::from_coefficients_vec(vec![
                -point.1.clone(),
                P::ScalarField::one(),
            ]);
        let mut quotient_polynomial_z1_y_coeffs = quotient_polynomial_y.coeffs.to_vec();
        quotient_polynomial_z1_y_coeffs.resize(y_srs.len(), <P::ScalarField>::zero());

        let proof = (
            P::G1::msm(&extended_powers, &extended_coeff).unwrap(), 
            P::G1::msm(&y_srs, &quotient_polynomial_z1_y_coeffs).unwrap());
        
        Ok(proof)
    }

    pub fn verify(
        v_srs: &VerifierSRS<P>,
        com: &P::G1,
        point: &(P::ScalarField, P::ScalarField),
        eval: &P::ScalarField,
        proof: &(P::G1, P::G1),
    ) -> Result<bool, Error> {
        let (x, y) = point;
        let left = P::pairing(com.clone() - v_srs.g * eval, v_srs.h.clone());
        let right1 = P::pairing(proof.0.clone(), v_srs.h_alpha.clone() - v_srs.h * x);
        let right2 = P::pairing(proof.1.clone(), v_srs.h_beta.clone() - v_srs.h * y);

        Ok(left == right1 + right2)
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use ark_bls12_381::Bls12_381;
    use ark_std::rand::{rngs::StdRng, SeedableRng};
    use ark_poly::DenseUVPolynomial;
    // use blake2::Blake2b;

    const BIVARIATE_X_DEGREE: usize = 8;
    const BIVARIATE_Y_DEGREE: usize = 6;
    //const UNIVARIATE_DEGREE: usize = 56;
    //const UNIVARIATE_DEGREE: usize = 65535;
    //const UNIVARIATE_DEGREE: usize = 1048575;

    type TestBivariatePolyCommitment = BivariateBatchKZG<Bls12_381>;
    // type TestUnivariatePolyCommitment = UnivariatePolynomialCommitment<Bls12_381, Blake2b>;

    #[test]
    fn bivariate_poly_commit_test() {
        let mut rng = StdRng::seed_from_u64(0u64);
        let srs =
            TestBivariatePolyCommitment::setup(&mut rng, BIVARIATE_X_DEGREE, BIVARIATE_Y_DEGREE)
                .unwrap();
        // let v_srs = srs.0.get_verifier_key();

        let mut x_polynomials = Vec::new();
        for _ in 0..BIVARIATE_Y_DEGREE + 1 {
            let mut x_polynomial_coeffs = vec![];
            for _ in 0..BIVARIATE_X_DEGREE + 1 {
                x_polynomial_coeffs.push(<Bls12_381 as Pairing>::ScalarField::rand(&mut rng));
            }
            x_polynomials.push(UnivariatePolynomial::from_coefficients_slice(
                &x_polynomial_coeffs,
            ));
        }
        let bivariate_polynomial = BivariatePolynomial { x_polynomials };

        // Commit to polynomial
        let com =
            TestBivariatePolyCommitment::commit(&srs.0, &bivariate_polynomial).unwrap();

        // Evaluate at challenge point
        let point = (UniformRand::rand(&mut rng), UniformRand::rand(&mut rng));
        let eval_proof = TestBivariatePolyCommitment::open(
            &srs.0,
            &bivariate_polynomial,
            &point
        )
        .unwrap();
        let eval = bivariate_polynomial.evaluate(&point);

        // proof size
        println!("Proof size is {} bytes", size_of_val(&eval_proof));

        // Verify proof
        assert!(
            TestBivariatePolyCommitment::verify(&srs.1, &com, &point, &eval, &eval_proof).unwrap()
        );

    }

}
