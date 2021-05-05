use num_traits::Float;

pub fn gamma<T: Float>(n: T) -> T {
    (n * T::epsilon()) / (T::one() - n * T::epsilon())
}
