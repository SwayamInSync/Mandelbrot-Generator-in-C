# Comparing Mandelbrot set for `quad-precision`, `double` and `long double` data types
Theoretical reference is taken from [Plotting algorithms for the Mandelbrot set(Wikipedia)](https://en.m.wikipedia.org/wiki/Plotting_algorithms_for_the_Mandelbrot_set)

## Updates
- [11/08/2024] Added Mandelbrot script in Python using Numpy with custom built QuadPrecision Datatype
- [11/08/2024] Added the `tlfloat` backend support for better quadprecision maths
  
## Comparison Table
| double | long double | __float128 |
|--------|-------------|------------|
| ![double](https://github.com/user-attachments/assets/4088a436-0057-4836-8b4d-31d9ef70a3d5) | ![long double](https://github.com/user-attachments/assets/738671ca-be91-4f24-a570-fad8307e6e6d) | ![__float128](https://github.com/user-attachments/assets/5f35b191-5339-445c-a19d-e2a4bf20b8b2) |


## Installation
- Compile and build [Sleef](https://sleef.org/)

## Usage
```
python setup.py build_ext --inplace
```
## Examples
> Mandelbrot Set
  
  ![normal](https://github.com/user-attachments/assets/868d6f4f-4dc6-422a-b83c-a25a9787ea54)
> SeaHorse Valley
  
  ![seahorse valley](https://github.com/user-attachments/assets/8642a51f-506c-46a4-b664-90e42454c0b3)
> SeaHorse Valley detailed
  
  ![seahorse valley detailed](https://github.com/user-attachments/assets/3fdf927f-a820-4de7-a0d6-fc2303c6c0cf)
> Spiral
  
  ![spiral](https://github.com/user-attachments/assets/4c7bb5e9-c45a-444e-bde1-c7b30fae45d4)
> Double Spiral
  
  ![double spiral](https://github.com/user-attachments/assets/f53541bc-0ac4-4d04-920f-6afb7552dbde)

> Spiral In The Bulb

  ![spiral in the bulb](https://github.com/user-attachments/assets/4f9eb262-4e91-41b5-b2ba-06679fa525c9)

> Elephant Valley

  ![elephant valley](https://github.com/user-attachments/assets/26cadd92-0ae8-47a1-9d0a-88371218084a)

> Embedded Julia Set

  ![embedded julia set](https://github.com/user-attachments/assets/f21e80ea-853a-4b31-852e-99476c04cb62)

> Fibonacci Spiral

  ![fibonacci spiral](https://github.com/user-attachments/assets/0e4f121d-c685-4292-85cb-a9cddf46728a)

> Intricate Filaments

  ![Intricate filaments](https://github.com/user-attachments/assets/851f252d-3aa3-4e3c-aa92-0d964dfff68e)

> Antenna

  ![Antenna](https://github.com/user-attachments/assets/0037eccd-b514-42bb-a81a-a2dca2dcd07c)



  



