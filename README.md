<a href="https://dearjohnsonny.github.io/Notes-from-DearJohn/">Notes from DearJohn</a>

<a href="https://dearjohnsonny.github.io/Notes-from-DearJohn-2/">Notes from DearJohn-2</a>

<a href="https://dearjohnsonny.github.io/Notes-from-DearJohn-3/">Notes from DearJohn-3</a>

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/195361749-01b1343d-0dc6-497b-bbc1-1b0ae149ca5e.png" width="900">
</div>

# Introduction to Vectors and Matrics

## Lengths and Dot Products
The length $\|v\|$ of a vector v is the square root of $\boldsymbol{v} \cdot \boldsymbol{v}$ 

$$
\text { length }=\|\boldsymbol{v}\|=\sqrt{\boldsymbol{v} \cdot \boldsymbol{v}}=\left(v_1^2+v_2^2+\cdots+v_n^2\right)^{1 / 2}
$$

$\boldsymbol{v} \cdot \boldsymbol{v}$ 相当于 $V^T V$

## Inverse Matrices

1 If the square matrix $A$ has an inverse, then both $A^{-1} A=I$ and $A A^{-1}=I$.

2 The algorithm to test invertibility is elimination: $A$ must have $n$ (nonzero) pivots.

3 The algebra test for invertibility is the determinant of $A$ : $\operatorname{det} A$ must not be zero.

4 The equation that tests for invertibility is $A \boldsymbol{x}=0: \boldsymbol{x}=\mathbf{0}$ must be the only solution.

5 If $A$ and $B$ (same size) are invertible then so is $A B:(A B)^{-1}=B^{-1} A^{-1}$.

6 $A A^{-1}=I$ is $n$ equations for $n$ columns of $A^{-1}$. Gauss-Jordan eliminates $[A I]$ to $\left[I A^{-1}\right]$.

7 The last page of the book gives 14 equivalent conditions for a square $A$ to be invertible.
