use rpbr::foundation::geometry::vector::Vector3;

fn main() {
    println!("Hello, world!");

    let vector = Vector3 {
        x: 1.0,
        y: 1.0,
        z: 2.0,
    };
    // println!("{}", vector);
    println!("{}", vector.length_squared());
}
