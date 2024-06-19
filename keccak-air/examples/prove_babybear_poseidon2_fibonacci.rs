use ark_ff::{BigInteger, PrimeField};
use p3_baby_bear::{BabyBear, DiffusionMatrixBabybear};
use p3_challenger::DuplexChallenger;
use p3_commit::ExtensionMmcs;
use p3_dft::Radix2DitParallel;
use p3_field::extension::BinomialExtensionField;
use p3_field::{AbstractField, Field};
use p3_fri::{FriConfig, TwoAdicFriPcs};
use p3_keccak_air::{FibonacciAir, FibonacciCols, TraceBorrowMut, NUM_FIBONACCI_COLS};
use p3_matrix::dense::RowMajorMatrix;
use p3_matrix::Matrix;
use p3_merkle_tree::FieldMerkleTreeMmcs;
use p3_poseidon2::Poseidon2;
use p3_symmetric::{PaddingFreeSponge, TruncatedPermutation};
use p3_uni_stark::{get_log_quotient_degree, prove, verify, StarkConfig, VerificationError};
use p3_util::log2_ceil_usize;
use rand::thread_rng;
use tracing_forest::util::LevelFilter;
use tracing_forest::ForestLayer;
use tracing_subscriber::layer::SubscriberExt;
use tracing_subscriber::util::SubscriberInitExt;
use tracing_subscriber::{EnvFilter, Registry};
use bincode;


fn main() -> Result<(), VerificationError> {
    let env_filter = EnvFilter::builder()
        .with_default_directive(LevelFilter::INFO.into())
        .from_env_lossy();

    Registry::default()
        .with(env_filter)
        .with(ForestLayer::default())
        .init();

    type Val = BabyBear;
    type Challenge = BinomialExtensionField<Val, 4>;

    type Perm = Poseidon2<Val, DiffusionMatrixBabybear, 16, 7>;
    

    let perm = Perm::new_from_rng(
        8,
        13,
        DiffusionMatrixBabybear,
        &mut thread_rng(),
    );

    type MyHash = PaddingFreeSponge<Perm, 16, 8, 8>;
    let hash = MyHash::new(perm.clone());

    type MyCompress = TruncatedPermutation<Perm, 2, 8, 16>;
    let compress = MyCompress::new(perm.clone());

    type ValMmcs = FieldMerkleTreeMmcs<
        <Val as Field>::Packing,
        <Val as Field>::Packing,
        MyHash,
        MyCompress,
        8,
    >;
    let val_mmcs = ValMmcs::new(hash, compress);

    type ChallengeMmcs = ExtensionMmcs<Val, Challenge, ValMmcs>;
    let challenge_mmcs = ChallengeMmcs::new(val_mmcs.clone());

    type Dft = Radix2DitParallel;
    let dft = Dft {};

    type Challenger = DuplexChallenger<Val, Perm, 16>;

    // 0..3
    // 3..6
    // 1 1 2
    // 1 2 3
    // ...
    const NUM_FIBONACCI_ROWS: usize = 1<<7;
    let mut trace =  RowMajorMatrix::new(vec![Val::zero(); NUM_FIBONACCI_ROWS * NUM_FIBONACCI_COLS], NUM_FIBONACCI_COLS);
    let rows = trace.borrow_rows_mut::<FibonacciCols<Val>>();
    rows[0].a = Val::from_canonical_u64(1);
    rows[0].b = Val::from_canonical_u64(1);
    rows[0].c = Val::from_canonical_u64(2);

    for i in 1..NUM_FIBONACCI_ROWS {
        rows[i].a = rows[i -1].b;
        rows[i].b = rows[i-1].c;
        rows[i].c = rows[i].a + rows[i].b ;
    }

    let fri_config = FriConfig {
        log_blowup: 1,
        num_queries: 100,
        proof_of_work_bits: 16,
        mmcs: challenge_mmcs,
    };
    type Pcs = TwoAdicFriPcs<Val, Dft, ValMmcs, ChallengeMmcs>;

    dbg!(log2_ceil_usize(trace.height()));
    dbg!(get_log_quotient_degree::<Val, FibonacciAir>(
        &FibonacciAir {},
        0
    ));

    let pcs = Pcs::new(log2_ceil_usize(trace.height()), dft, val_mmcs, fri_config);

    type MyConfig = StarkConfig<Pcs, Challenge, Challenger>;
    let config = MyConfig::new(pcs);

    let mut challenger = Challenger::new(perm.clone());

    let proof = prove::<MyConfig, _>(&config, &FibonacciAir {}, &mut challenger, trace, &vec![]);

    println!("proof size: {} bytes", bincode::serialize(&proof).unwrap().len());

    std::fs::write(
        "proof_poseidon2_fibonacci.json",
        serde_json::to_string(&proof).unwrap(),
    )
    .unwrap();

    let mut challenger = Challenger::new(perm);
    verify(&config, &FibonacciAir {}, &mut challenger, &proof, &vec![]).unwrap();
    Ok(())
}
