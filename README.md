<a href="https://dearjohnsonny.github.io/Notes1-Biotech/">Notes1-Biotech</a>

<a href="https://dearjohnsonny.github.io/Notes2-Biotech/">Notes2-Biotech</a>

<a href="https://dearjohnsonny.github.io/Notes3-Statistics/">Notes3-Statistics</a>

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/196834390-53f62718-18ed-4975-b607-bcaa4911e254.png" width="1500">
</div>


<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/196833895-0f60f1ee-8163-4ae3-8eef-8c509937e251.png" width="1500">
</div>

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/196833845-03c478c0-177f-40a0-95ab-f19df39c91d1.png" width="1500">
</div>


# Determinants 行列式
* 方阵才有行列式，行列式的结果是一个值（可正可负）。
* 行列式为0的矩阵为奇异矩阵，不可逆（没有乘法逆元）
* 当某个方阵的行列式结果为0时，证明该方阵表示的线性变换将把原基向量压缩到更低维的空间（低一维还是二维不清楚）。因为该线性变换的（空间的）基是线性相关的，因此某些基可以用其他基表示，因此该线性变换至少损失一个基底，维度至少降低一维。
  ![image](https://user-images.githubusercontent.com/111955215/195834172-03c9b353-145b-4824-bdd3-4706b9b2651c.png)
  如上图中的两个新基底，实际上是共线的，也就是线性相关，因此必定会损失y轴方向的基地，从2维变成了1维
* 行列式的值的绝对值表示该（改变原基向量）的线性变换形成的新基底的乘积放大（或缩小）多少倍（如果是2 × 2矩阵，则是拉伸x轴和y轴的基底并产生的新四边形的面积；若为3 × 3矩阵，则为体积）
  ![image](https://user-images.githubusercontent.com/111955215/195833908-c2efe3d4-5124-4a43-9a12-4f85aec2e669.png)
* 该结果的正负表示空间上的变向（如二维的翻面）
* 借助下图理解一下行列式的展开
  ![image](https://user-images.githubusercontent.com/111955215/195835018-ada169e3-a73d-4afc-a5bd-760b56830d0e.png)

## 克莱姆法则
在引入克莱姆法则之前，先引入有关n元线性方程组和有关矩阵、行列式的概念。含有n个未知数的线性方程组称为n元线性方程组。

$$
\left\{\begin{array}{l}
a_{11} x_1+a_{12} x_2+\cdots+a_{1 n} x_n=b_1 \\
a_{21} x_1+a_{22} x_2+\cdots+a_{2 n} x_n=b_2 \\
\cdots \cdots \cdots \cdots \cdots \cdots \cdots \cdots \cdots \\
a_{n 1} x_1+a_{n 2} x_2+\cdots+a_{n n} x_n=b_n
\end{array}\right.
$$

当常数项全为零时，线性方程组⑵称为齐次线性方程组：

$$
\left\{\begin{array}{l}
a_{11} x_1+a_{12} x_2+\cdots+a_{1 n} x_n=0 \\
a_{21} x_1+a_{22} x_2+\cdots+a_{2 n} x_n=0 \\
\cdots \cdots \cdots \cdots \cdots \cdots \cdots \cdots \\
a_{n 1} x_1+a_{n 2} x_2+\cdots+a_{n n} x_n=0
\end{array}\right.
$$

系数构成的行列式称为该方程组的系数行列式D：

$$
\text { 即 } D=\left|\begin{array}{cccc}
a_{11} & a_{12} & \cdots & a_{1 n} \\
a_{21} & a_{22} & \cdots & a_{2 n} \\
\vdots & \vdots & \ddots & \vdots \\
a_{n 1} & a_{n 2} & \cdots & a_{n n}
\end{array}\right|
$$

定理
记法1: 若线性方程组(1)的系数矩阵可逆 (非奇异)，即系数行列式 $\mathrm{D} \neq 0$ 。有唯一解，其解为 $X_0=A^{-1} \beta$

记法2: 若线性方程组 (1)的系数矩阵可逆（非奇异），即系数行列式 $D \neq 0$ ，则线性方程组(1)有唯一解，其解为 $x_j=\frac{D_j}{D}(j=1,2, \cdots, n)$

**其中 $D_j$ 是把D中第j列元素对应地换成常数项而其余各列保持不变所得到的行列式。**

记法1是将解写成矩阵 (列向量) 形式，而记法 2 是将解分别写成数字，本质相同。当其右端的常数项b1,b2,...,bn不全为零时，线性方程组⑴称为非齐次线性方程组。

## 空间几何

向量的外积：$\left.\left(x_1, y_1, z_1\right) \times\left(x_2, y_2, z_2\right)=\left(\left|\begin{array}{ll}y_1 & z_1 \\ y_2 & z_2\end{array}\right|, \mid \begin{array}{cc}z_1 & x_1 \\ z_2 & x_2\end{array}\right],\left|\begin{array}{ll}x_1 & y_1 \\ x_2 & y_2\end{array}\right|\right)$


# Introduction to Vectors and Matrics

## Lengths and Dot Products
The length $\|v\|$ of a vector v is the square root of $\boldsymbol{v} \cdot \boldsymbol{v}$ 

$$
\text { length }=\|\boldsymbol{v}\|=\sqrt{\boldsymbol{v} \cdot \boldsymbol{v}}=\left(v_1^2+v_2^2+\cdots+v_n^2\right)^{1 / 2}
$$

$\boldsymbol{v} \cdot \boldsymbol{v}$ 相当于 $V^T V$

## Inverse Matrices

1 If the square matrix $A$ has an inverse, then both $A^{-1} A=I$ and $A A^{-1}=I$

2 The algorithm to test invertibility is elimination: $A$ must have $n$ (nonzero) pivots.

3 The algebra test for invertibility is the determinant of $A$ : $\operatorname{det} A$ must not be zero.

4 The equation that tests for invertibility is $A \boldsymbol{x}=0: \boldsymbol{x}=\mathbf{0}$ must be the only solution.

5 If $A$ and $B$ (same size) are invertible then so is $A B:(A B)^{-1}=B^{-1} A^{-1}$.

6 $A A^{-1}=I$ is $n$ equations for $n$ columns of $A^{-1}$. Gauss-Jordan eliminates $[A I]$ to $\left[I A^{-1}\right]$.

7 The last page of the book gives 14 equivalent conditions for a square $A$ to be invertible.

## A=LU分解法

对2×2矩阵，要得到U这个上三角矩阵，需要我们将2行1列的元消去，则需要E21这个矩阵，然后再求逆矩阵即可

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/195482901-fb755dad-6ca3-433e-a835-85b4726ac0de.png" width="500">
</div>

对于更广泛的情况（如三维方阵），也是一样的思路：

$$
\left(E_{32} E_{31} E_{21}\right) A=U \text { becomes } A=\left(E_{21}^{-1} E_{31}^{-1} E_{32}^{-1}\right) U \text { which is } A=L U \text {. }
$$

**Better balance from LDU**：A = L U is "unsymmetric" because U has the pivots on its diagonal where L has l's. This is easy to change. Divide U by a diagonal matrix D that contains the pivots. That leaves a new triangular matrix with l's on the diagonal(对角线): 

$$
\text { The triangular factorization can be written } A=L U \text { or } A=L D U \text {. }
$$

如下图的例子：

$$
\left[\begin{array}{ll}
1 & 0 \\
3 & 1
\end{array}\right]\left[\begin{array}{ll}
2 & 8 \\
0 & 5
\end{array}\right] \quad \text { splits further into } \quad\left[\begin{array}{ll}
1 & 0 \\
3 & 1
\end{array}\right]\left[\begin{array}{ll}
2 & \\
& 5
\end{array}\right]\left[\begin{array}{ll}
1 & 4 \\
0 & 1
\end{array}\right]
$$

## Transposes and Permutations转置与排列
1 The transposes of $A \boldsymbol{x}$ and $A B$ and $A^{-1}$ are $\boldsymbol{x}^{\mathrm{T}} A^{\mathrm{T}}$ and $B^{\mathrm{T}} A^{\mathrm{T}}$ and $\left(A^{\mathrm{T}}\right)^{-1}$.

2 The dot product (inner product) is $\boldsymbol{x} \cdot \boldsymbol{y}=\boldsymbol{x}^{\mathrm{T}} \boldsymbol{y}$. This is $(1 \times n)(n \times 1)=(1 \times 1)$. The outer product is $\boldsymbol{x} \boldsymbol{y}^{\mathbf{T}}=$ column times row $=(n \times 1)(1 \times n)=n \times n$ matrix .

3 The idea behind $A^{\mathrm{T}}$ is that $A \boldsymbol{x} \cdot \boldsymbol{y}$ equals $\boldsymbol{x} \cdot A^{\mathrm{T}} \boldsymbol{y}$ because $(A \boldsymbol{x})^{\mathrm{T}} \boldsymbol{y}=\boldsymbol{x}^{\mathrm{T}} A^{\mathrm{T}} \boldsymbol{y}=\boldsymbol{x}^{\mathrm{T}}\left(A^{\mathrm{T}} \boldsymbol{y}\right)$.

4 A symmetric matrix has $S^{\mathrm{T}}=\boldsymbol{S}$ (and the product $A^{\mathrm{T}} A$ is always symmetric).

5 An orthogonal matrix has $\boldsymbol{Q}^{\mathrm{T}}=\boldsymbol{Q}^{-1}$. The columns of $Q$ are orthogonal unit vectors.

6 A permutation matrix $\boldsymbol{P}$ has the same rows as $I$ (in any order). There are $\boldsymbol{n}$ ! different orders.

7 Then $P x$ puts the components $x_1, x_2, \ldots, x_n$ in that new order. And $P^{\mathrm{T}}$ equals $P^{-1}$.

$$
\text { The transpose of } A^{\mathrm{T}} A \text { is } A^{\mathrm{T}}\left(A^{\mathrm{T}}\right)^{\mathrm{T}} \text { which is } A^{\mathrm{T}} A \text { again. }
$$

The matrix $A A^{\mathrm{T}}$ is also symmetric，for instance: 

$$
A A^{\mathrm{T}}=\left[\begin{array}{rr}
2 & -1 \\
-1 & 2
\end{array}\right] \text { and } A^{\mathrm{T}} A=\left[\begin{array}{rrr}
1 & -1 & 0 \\
-1 & 2 & -1 \\
0 & -1 & 1
\end{array}\right]
$$

If $S=S^{\mathrm{T}}$ is factored into $L D U$ with no row exchanges, then $U$ is exactly $L^{\mathrm{T}}$.

The symmetric factorization of a symmetric matrix is $S=L D L^{\mathrm{T}}$.

# Spaces of Vectors 

M The vector space of all real 2 by 2 matrices.

F The vector space of all real functions $f(x)$.

$\mathrm{Z}$ The vector space that consists only of a zero vector.

如何在向量空间中选择基底：

The column space of $A$-choose the pivot columns of $A$ as a basis.

The row space of $A$-choose the nonzero rows of $R$ as a basis.

The nullspace of $A-$ choose the special solutions to $R \boldsymbol{x}=\mathbf{0}$ (and $A \boldsymbol{x}=\mathbf{0}$ ).

## 关于列空间column space
**The column space** consists of all linear combinations of the columns. The combinations are all possible vectors $A x$. They fill the column space $C(A)$.

The system $A x=b$ is solvable if and only if $b$ is in the column space of $A$.

列空间是一种用向量空间中的向量来张成子空间的方法，但可以将其推广：
**Important** Instead of columns in $\mathbf{R}^m$, we could start with any set $\mathbf{S}$ of vectors in a vector space $\mathbf{V}$. To get a subspace $\mathbf{S S}$ of $\mathbf{V}$, we take all combinations of the vectors in that set:

$$
\begin{gathered}
\mathbf{S}=\text { set of vectors in } \mathbf{V} \text { (probably not a subspace) } \\
\mathbf{S S}=\text { all combinations of vectors in } \mathbf{S} \text { (definitely a subspace) } \\
\mathbf{S S}=\text { all } c_1 \boldsymbol{v}_1+\cdots+c_N \boldsymbol{v}_N=\text { the subspace of } \mathbf{V} \text { "spanned" by } \mathbf{S}
\end{gathered}
$$

## 关于零空间The Nullspace of A: Solving Ax = 0 and Rx = 0
1 The nullspace $\boldsymbol{N}(A)$ in $\mathbf{R}^n$ contains all solutions $\boldsymbol{x}$ to $A \boldsymbol{x}=\mathbf{0}$. This includes $\boldsymbol{x}=\mathbf{0}$.

2 Elimination (from $A$ to $U$ to $R$ ) does not change the nullspace: $\boldsymbol{N}(A)=\boldsymbol{N}(U)=\boldsymbol{N}(R)$.

3 The reduced row echelon form $R=\operatorname{rref}(A)$ has all pivots $=1$, with zeros above and below.

4 If column $j$ of $R$ is free (no pivot), there is a "special solution" to $A \boldsymbol{x}=\mathbf{0}$ with $x_j=1$.

5 Number of pivots $=$ number of nonzero rows in $R=\operatorname{rank} \boldsymbol{r}$. There are $n-r$ free columns.

6 Every matrix with $m < n$ has nonzero solutions to $A \boldsymbol{x}=\mathbf{0}$ in its nullspace.

The rank of $A$ is the number of pivots. This number is $r$ 下图矩阵的rank为3

$$
\boldsymbol{R}=\left[\begin{array}{lllllll}
\mathbf{1} & \mathbf{0} & x & x & x & \mathbf{0} & x \\
0 & 1 & x & x & x & 0 & x \\
\mathbf{0} & \mathbf{0} & 0 & 0 & 0 & \mathbf{1} & x \\
0 & 0 & 0 & 0 & 0 & 0 & 0
\end{array}\right] \quad \begin{aligned}
&\text { Three pivot variables } x_1, x_2, x_6 \\
&\text { Four free variables } x_3, x_4, x_5, x_7 \\
&\text { Four special solutions } s \text { in } N(\boldsymbol{R}) \\
&\text { The pivot rows and columns contain } \boldsymbol{}
\end{aligned}
$$

## The Complete Solution to Ax = b

1 Complete solution to $A \boldsymbol{x}=\boldsymbol{b}: \boldsymbol{x}=$ (one particular solution $\left.\boldsymbol{x}_p\right)+\left(\right.$ any $\boldsymbol{x}_n$ in the nullspace).

2 Elimination on $\left[\begin{array}{ll}A & \boldsymbol{b}\end{array}\right]$ leads to $\left[\begin{array}{ll}R & \boldsymbol{d}\end{array}\right]$. Then $A \boldsymbol{x}=\boldsymbol{b}$ is equivalent to $R \boldsymbol{x}=\boldsymbol{d}$.

3 A \boldsymbol{x}=\boldsymbol{b}$ and $R \boldsymbol{x}=\boldsymbol{d}$ are solvable only when all zero rows of $R$ have zeros in $\boldsymbol{d}$.

4 When $R \boldsymbol{x}=\boldsymbol{d}$ is solvable, one very particular solution $\boldsymbol{x}_p$ has all free variables equal to zero.

5 has full column rank $\boldsymbol{r}=\boldsymbol{n}$ when its nullspace $\boldsymbol{N}(A)=$ zero vector: no free variables.

6 A has full row rank $\boldsymbol{r}=\boldsymbol{m}$ when its column space $\boldsymbol{C}(A)$ is $\mathbf{R}^m: A \boldsymbol{x}=\boldsymbol{b}$ is always solvable.

7 The four cases are $r=m=n$ ( $A$ is invertible) and $r=m < n$ (every $A \boldsymbol{x}=\boldsymbol{b}$ is solvable) and $r=n < m(A \boldsymbol{x}=\boldsymbol{b}$ has 1 or 0 solutions) and $r < m, r < n$ ( 0 or $\infty$ solutions).

8 $A \boldsymbol{x}=\boldsymbol{b}$ is solvable if and only if the last $m-r$ equations reduce to $0=0$

要解Ax = b这个方程，需要一个AX = b方程的特解，还需要一个AX = 0的方程的通解：

$$
\begin{array}{lll}
x_{\text {particular }} & \text { The particular solution solves } & A x_p=b \\
x_{\text {nullspace }} & \text { The } n-r \text { special solutions solve } & A x_n=0 .
\end{array}
$$

比如下面这个例子：

$$
R x_p=\left[\begin{array}{llll}
1 & 3 & \mathbf{0} & 2 \\
\mathbf{0} & 0 & \mathbf{1} & 4 \\
0 & 0 & 0 & 0
\end{array}\right]\left[\begin{array}{l}
\mathbf{1} \\
0 \\
\mathbf{6} \\
0
\end{array}\right]=\left[\begin{array}{ll}
\mathbf{1} \\
\mathbf{6} \\
0
\end{array}\right] \quad \begin{aligned}
&\text { Pivot variables } \mathbf{1}, \mathbf{6} \\
&\text { Solution } x_p=(\mathbf{1}, \mathbf{0}, \mathbf{6}, \mathbf{0})
\end{aligned}
$$

解特解时，将自由变量设定为0，则其中一个特解肯定是[1,0,6,0]^T，AX = 0的通解分别令自由变量为0和1，解出来之后再乘上相应的Xi即可：

$$
\boldsymbol{x}=x_p+\boldsymbol{x}_n=\left[\begin{array}{l}
1 \\
0 \\
6 \\
0
\end{array}\right]+x_2\left[\begin{array}{r}
-3 \\
1 \\
0 \\
0
\end{array}\right]+x_4\left[\begin{array}{r}
-2 \\
0 \\
-4 \\
1
\end{array}\right]
$$

再如下例：

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/195583662-ba35f396-18a1-40b4-a6b2-5fa30031b9d2.png" width="800">
</div>

### 列满秩与行满秩
* The rank r = n，则该矩阵A应该是瘦长型的(m ≥ n)，将A简化成为R将得到下式：

$$
R=\left[\begin{array}{l}
I \\
0
\end{array}\right]=\left[\begin{array}{l}
n \text { by } n \text { identity matrix } \\
m-n \text { rows of zeros }
\end{array}\right]
$$

Every matrix $A$ with full column $\operatorname{rank}(\boldsymbol{r}=\boldsymbol{n})$ has all these properties:
1. All columns of $A$ are pivot columns.
2. There are no free variables or special solutions.
3. The nullspace $\boldsymbol{N}(A)$ contains only the zero vector $\boldsymbol{x}=\mathbf{0}$.
4. If $A \boldsymbol{x}=b$ has a solution (it might not) then it has only one solution.

* A matrix has full row rank if r = m，则该矩阵是矮胖型的(n ≥ m)
Every matrix $A$ with full row $\operatorname{rank}(\boldsymbol{r}=\boldsymbol{m})$ has all these properties:
1. All rows have pivots, and $R$ has no zero rows.
2. $A \boldsymbol{x}=\boldsymbol{b}$ has a solution for every right side $b$.
3. The column space is the whole space $\mathbf{R}^m$.
4. There are $n-r=n-m$ special solutions in the nullspace of $A$.

要解方程则按下图：

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/195591903-ce4c2bbd-199b-4d1e-b571-7a92b1bd03e8.png" width="800">
</div>

### 几个基本子空间
#### 首先N(A)⊥R(A^T)
令 $A$ 为一 $m \times n$ 矩阵, 并令 $\boldsymbol{x} \in N(A), N(A)$ 为 $A$ 的零空间. 由于 $A \boldsymbol{x}=\mathbf{0}$, 我们有 $a_{i 1} x_1+a_{i 2} x_2+\cdots+a_{i n} x_n=0$
其中 $i=1, \cdots, m$. 方程 (1) 说明, $\boldsymbol{x}$ 与 $A^{\mathrm{T}}$ 的第 $i$ 个列向量正交, 其中 $i=1, \cdots, m$. 由 于 $\boldsymbol{x}$ 和 $A^{\mathrm{T}}$ 的每一个列向量正交, 所以它和 $A^{\mathrm{T}}$ 的列向量的任何线性组合也正交. 因此, 若 $\boldsymbol{y}$ 为 $A^T$ 的列空间中的任何一个向量, 则 $\boldsymbol{x}^{\mathrm{T}} \boldsymbol{y}=0$. 于是, $N(A)$ 中的每一向量都和 $A^{\mathrm{T}}$ 的列空间中的任何向量正交. 当 $\mathbf{R}^n$ 的两个子空间具有这个性质( $\boldsymbol{x}^{\mathrm{T}} \boldsymbol{y}=0$ )时, 称它们是正交的.

#### 正交补
![image](https://user-images.githubusercontent.com/111955215/195830293-ff4aabae-5d61-4f8e-934b-b001d6f4fff7.png)

正交补的概念要考虑**取满**能正交的集合

#### N(A)与R(A^T)互为正交补（基本子空间定理）
基本子空间定理：若 $A$ 为 $-m \times n$ 矩阵, 则 $N(A)=R\left(A^{\mathrm{T}}\right)^{\perp}$, 且 $N\left(A^{\mathrm{T}}\right)=R(A)^{\perp}$.

#### 互为正交补的两个空间维数互补

想象一下三维空间中的一维和二维空间则很好理解这个定理

# 正交性 Orthogonality 
## 线性无关、基与维数
1. Independent vectors(no extra vectors) 
2. Spanning a space(enough vectors to produce the rest) 
3. Basis for a space(not too many or too few) 少了不够，多了不独立
4. Dimension of a space(the number of vectors in a basis) 


1 Independent columns of $A$ : The only solution to $A \boldsymbol{x}=\mathbf{0}$ is $\boldsymbol{x}=\mathbf{0}$. The nullspace is $\boldsymbol{Z}$.

2 Independent vectors: The only zero combination $c_1 \boldsymbol{v}_1+\cdots+c_k \boldsymbol{v}_k=0$ has all $c$ 's $=0$.

3 A matrix with $m < n$ has dependent columns : At least $n-m$ free variables/ special solutions.

4 The vectors $\boldsymbol{v}_1, \ldots, \boldsymbol{v}_k$ span the space $\boldsymbol{S}$ if $\boldsymbol{S}=$ all combinations of the $\boldsymbol{v}$ 's.

5 The vectors $\boldsymbol{v}_1, \ldots, \boldsymbol{v}_k$ are a basis for $\boldsymbol{S}$ if they are independent and they span $\boldsymbol{S}$.

6 The dimension of a space $S$ is the number of vectors in every basis for $S$.

7 If $A$ is 4 by 4 and invertible, its columns are a basis for $\mathbf{R}^4$. The dimension of $\mathbf{R}^4$ is 4 .

关于线性无关：当 $A \boldsymbol{x}=\mathbf{0}$ 的解也就是 $x_1 \boldsymbol{v}_1+x_2 \boldsymbol{v}_2+\cdots+x_n \boldsymbol{v}_n=\mathbf{0}$ 的解只有 $\boldsymbol{x}=0$时，可以确定A的列向量是线性无关的 

The columns are certainly dependent if n > m, because Ax = 0 has a nonzero solution（矮胖型的零空间必定不止包含原点；未知数多方程少，肯定可以解出来）

### 行空间
The row space of $A$ is $\boldsymbol{C}\left(A^{\mathrm{T}}\right)$. It is the column space of $A^{\mathrm{T}}$.This row space of A is a subspace of $\mathbf{R}^2$

也可以记做R(A)

## 投影projection

When $b$ is projected onto a line, its projection $p$ is the part of $b$ along that line. If $b$ is projected onto a plane, $p$ is the part in that plane. The projection $p$ is $P b$.其中P为投影矩阵

##$When $b$ is projected onto a line, its projection $p$ is the part of $b$ along that line. If $b$ is projected onto a plane, $p$ is the part in that plane. The projection $p$ is $P b$.

### 投影到线上
![image](https://user-images.githubusercontent.com/111955215/195606114-f9c3a016-7e75-475e-ae3e-fa107101f575.png)

以上图为例：b为原始向量，a为目标线，p为投影向量，P为投影矩阵

Projecting $b$ onto $\boldsymbol{a}$ with error $\boldsymbol{e}=\boldsymbol{b}-\widehat{x} \boldsymbol{a}$ $\boldsymbol{a} \cdot(\boldsymbol{b}-\widehat{\boldsymbol{x}} \boldsymbol{a})=0 \quad$ or $\quad \boldsymbol{a} \cdot \boldsymbol{b}-\widehat{\boldsymbol{x}} \boldsymbol{a} \cdot \boldsymbol{a}=0$

$$
\widehat{\boldsymbol{x}}=\frac{\boldsymbol{a} \cdot \boldsymbol{b}}{\boldsymbol{a} \cdot \boldsymbol{a}}=\frac{\boldsymbol{a}^{\mathrm{T}} \boldsymbol{b}}{\boldsymbol{a}^{\mathrm{T}} \boldsymbol{a}}
$$

Projection $\boldsymbol{p}=\boldsymbol{a} \widehat{\boldsymbol{x}}=\boldsymbol{a} \frac{\boldsymbol{a}^{\mathrm{T}} \boldsymbol{b}}{\boldsymbol{a}^{\mathrm{T}} \boldsymbol{a}}=P \boldsymbol{b}$ when the matrix is $\quad P=\frac{\boldsymbol{a a}^{\mathrm{T}}}{\boldsymbol{a}^{\mathrm{T}} \boldsymbol{a}}$.

### 投影到面上

$$
\begin{array}{cc}
\boldsymbol{a}_1^{\mathrm{T}}(\boldsymbol{b}-A \widehat{\boldsymbol{x}})=0 \\
\vdots \\
\boldsymbol{a}_n^{\mathrm{T}}(\boldsymbol{b}-A \widehat{\boldsymbol{x}})=0
\end{array} \quad \text { or } \quad\left[\begin{array}{c}
-\boldsymbol{a}_1^{\mathrm{T}}- \\
\vdots \\
-\boldsymbol{a}_n^{\mathrm{T}}-
\end{array}\right]\left[\begin{array}{l}
\boldsymbol{b}-A \widehat{\boldsymbol{x}}
\end{array}\right]=\left[\begin{array}{l}
\mathbf{0}
\end{array}\right]
$$

$$
\text { The matrix with those rows } \boldsymbol{a}_i^{\mathrm{T}} \text { is } A^{\mathrm{T}} \text {. The } n \text { equations are exactly } A^{\mathrm{T}}(\boldsymbol{b}-A \widehat{\boldsymbol{x}})=\mathbf{0} \text {. }
$$

The combination $\boldsymbol{p}=\widehat{x}_1 \boldsymbol{a}_1+\cdots+\widehat{x}_n \boldsymbol{a}_n=A \widehat{\boldsymbol{x}}$ that is closest to $\boldsymbol{b}$ comes from $\widehat{\boldsymbol{x}}$ : Find $\widehat{x}(n \times 1) \quad A^{\mathrm{T}}(\boldsymbol{b}-A \widehat{\boldsymbol{x}})=\mathbf{0} \quad$ or $\quad A^{\mathrm{T}} A \widehat{\boldsymbol{x}}=A^{\mathrm{T}} \boldsymbol{b}$
This symmetric matrix $A^{\mathrm{T}} A$ is $n$ by $n$. It is invertible if the $\boldsymbol{a}$ 's are independent. The solution is $\widehat{\boldsymbol{x}}=\left(A^{\mathrm{T}} A\right)^{-1} A^{\mathrm{T}} \boldsymbol{b}$. The projection of $\boldsymbol{b}$ onto the subspace is $\boldsymbol{p}$ :
Find $\boldsymbol{p}(m \times 1)$

$$
\boldsymbol{p}=A \widehat{\boldsymbol{x}}=A\left(A^{\mathrm{T}} A\right)^{-1} A^{\mathrm{T}} \boldsymbol{b} .
$$

The next formula picks out the projection matrix that is multiplying $b$ in:
Find $\boldsymbol{P}(m \times m)$

$$
P=A\left(A^{\mathrm{T}} A\right)^{-1} A^{\mathrm{T}} .
$$

# 特征值和特征向量
对一个空间进行线性变换之后发现大多数向量都会离开原来所张成的空间
![image](https://user-images.githubusercontent.com/111955215/195835359-26314d04-824f-4fe9-b4ea-b670585a5b25.png)

而有的特殊向量会留在原空间中，其被称为特征向量
![image](https://user-images.githubusercontent.com/111955215/195835643-9629d674-ddda-41f8-a88c-8686aa507114.png)

还比如三维空间中的旋转轴
![image](https://user-images.githubusercontent.com/111955215/195835891-a6600ba5-293f-425b-a305-a39ca591c959.png)

而该拉伸的倍数被称为特征值（该例下特征值为1）

想象特征值和特征向量的拉伸非常有助于我们理解该线性变换（如下图）
![image](https://user-images.githubusercontent.com/111955215/195836489-f71612a8-f244-4ab6-aa15-2127184a1d95.png)

要解特征向量就非常简单了，只需要将这个式子 $A \overrightarrow{\mathbf{v}}=\lambda \overrightarrow{\mathbf{v}}$ 做一下变换，向量乘对应的I为本身，从而将其变为两个矩阵相等 $A \overrightarrow{\mathbf{v}}=(\lambda I) \overrightarrow{\mathbf{v}}$ ，从而得到 $(A-\lambda I) \overrightarrow{\mathbf{v}}=\overrightarrow{0}$

而要一个非零向量进行左边的线性变换 $A-\lambda I$ 之后变为0向量，则其维数降低，也就是需要该线性变换将其压缩，因此 $\operatorname{det}(A-\lambda I)=0$

## 对角化 Diagonalizing 也叫作相似对角化

**对角化的结果是一个对角矩阵，本质就是把矩阵列向量都放到标准轴上。可以说，对角矩阵一定是“观看演出时”的最佳视角**

![image](https://user-images.githubusercontent.com/111955215/195865438-f35858ed-ba9d-4a94-af43-b7fa49afb657.png)

1 The columns of $A X=X \Lambda$ are $A \boldsymbol{x}_k=\lambda_k \boldsymbol{x}_k$. The eigenvalue matrix $\Lambda$ is diagonal.

2 $\boldsymbol{n}$ independent eigenvectors in $X$ diagonalize $A \quad A=\boldsymbol{X} \boldsymbol{\Lambda} \boldsymbol{X}^{-1}$ and $\boldsymbol{\Lambda}=\boldsymbol{X}^{-1} \boldsymbol{A} \boldsymbol{X}$

3 The eigenvector matrix $X$ also diagonalizes all powers $A^k$ : $A^k=\boldsymbol{X} \Lambda^k \boldsymbol{X}^{-1}$

4 No equal eigenvalues $\Rightarrow X$ is invertible and $A$ can be diagonalized.
Equal eigenvalues $\Rightarrow A$ might have too few independent eigenvectors. Then $X^{-1}$ fails.

5 Every matrix $C=B^{-1} A B$ has the same eigenvalues as $A$. These $C$ 's are "similar" to $A$.

### 正交矩阵必定可以对角化
正交矩阵的定义：称 $\mathrm{n}$ 阶方阵 $\mathrm{A}$ 是正交矩阵，若 $A^T A=I$

正交矩阵有几个重要性质:

1. A的逆等于 $A$ 的转置，即 $A^{-1}=A^T$ （显然）

2. $\mathrm{A}$ 的行列式为 $\pm 1$ ， 即 $|A|=\pm 1$ （通过二维的面积和三维的体积来理解这个事）

3. A的行 (列) 向量组为 $n$ 维单位正交向量组（矩阵的各列向量都是单位向量，并且两两正交；对于正交矩阵，组成它的列向量构成了一个空间的基，称之为：规范正交基。）

正交矩阵的各个基底本身就是相互垂直，只是说它不见得是各个标准轴，因为可能并不放在各个标准轴上，可能有移动或者旋转

**对角化的结果是一个对角矩阵，本质就是把矩阵列向量都放到标准轴上**。 那么很显然：正交矩阵一定可以做到！

因此：凡是正交矩阵一定可以对角化

## 需要提到的基变换
常规选择的基底是i和j，而有时会选择任意其他的b1和b2作为基底。如果希望对该b1和b2作为基底的某个向量进行线性变换，需要先将该基底变换为常规基底i和j，因为如果将矩阵视为线性变换的前提是基底为单位向量。

用i和j来理解b1和b2的话，只需要将两个基底分别作为列向量。下图即为将在b1和b2下表示为[-1,2]的向量转换为常规基底下的向量的过程：

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/195962958-d7d02cd3-746a-4bbc-a653-5df00da425cf.png" width="600">
</div>

因此可以得到下列的基变换过程，先用常规基底来刻画向量，再线性变换，最后再用逆矩阵将该变换还原到原来的基底：

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/195963110-c312fe9a-b705-47a2-9f22-f3b52f79e5c7.png" width="900">
</div>


需要提到的是：表达式 $A^{-1} M A$ 暗示着一种数学上的转移作用，中间的矩阵代表一种线性变换，而外侧的两个矩阵代表着转移作用

## 相似矩阵 Similar Matrices: Same Eigenvalues 

$$
\text { All the matrices } A=B C B^{-1} \text { are "similar." They all share the eigenvalues of } C \text {. }
$$

Proof 

Suppose $C \boldsymbol{x}=\lambda \boldsymbol{x}$. Then $B C B^{-1}$ has the same eigenvalue $\lambda$ with the new eigenvector $B \boldsymbol{x}$ :
Same $\boldsymbol{\lambda} \quad\left(B C B^{-1}\right)(B \boldsymbol{x})=B C \boldsymbol{x}=B \lambda \boldsymbol{x}=\lambda(B \boldsymbol{x})$

在选择了两个不同的角度（线性空间的两组基）之后，同一个线性变换所对应的矩阵是不一样的。而这些矩阵都是彼此相似的矩阵。也就是说，相似矩阵是同一个线性变换在不同基下的矩阵，因此线性变换相同，则对应的线性变换的特征值也相同，只是选择的基底不同而已。

只要找到这两个组基的过渡矩阵（实际上就是基变换矩阵），就可以轻松得到相似定义中的公式

### 相似矩阵具有相同的特征多项式 Similar matrices have the same characteristic polynomial.

PS：特征多项式是行列式的展开： $\operatorname{det}(x I-A)$ 
Proof. If $B=P^{-1} A P$, then

$$
\begin{aligned}
\operatorname{det}(x I-B) &=\operatorname{det}\left(x I-P^{-1} A P\right) \\
&=\operatorname{det}\left(P^{-1}(x I-A) P\right) \\
&=\operatorname{det} P^{-1} \cdot \operatorname{det}(x I-A) \cdot \operatorname{det} P \\
&=\operatorname{det}(x I-A)
\end{aligned}
$$

## 奇异值分解 Singular value decomposition

<a href="https://blog.csdn.net/nstarLDS/article/details/106206805">不错的证明</a>

<a href="https://zhuanlan.zhihu.com/p/29846048">可以看看的例子</a>

以2 × 2矩阵为例介绍奇异值分解的意义：

**Top**: The action of M, indicated by its effect on the unit disc D and the two canonical unit vectors e1 and e2.

**Left**: The action of V⁎, a rotation, on D, e1, and e2.

**Bottom**: The action of Σ, a scaling by the singular values σ1 horizontally and σ2 vertically.

**Right**: The action of U, another rotation.

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/195971385-543fc40b-410d-48c6-8933-f8c4e1bb8fe6.png" width="500">
</div>

从线性变换的角度理解奇异值分解, $m \times n$ 矩阵 $A$ 表示从 $n$ 维空间 $\mathbf{R}^n$ 到 $m$ 维空间 $\mathbf{R}^m$ 的一个线性变换，

$$
T: x \rightarrow A x
$$

$x \in \mathbf{R}^n, A x \in \mathbf{R}^m, x$ 和 $A x$ 分别是各自空间的向量。线性变换可以分解为三个简单 的变换: 一个坐标系的旋转或反射变换、一个坐标轴的缩放变换、另一个坐标系的旋 转或反射变换。奇异值定理保证这种分解一定存在。这就是奇异值分解的几何解释。

对矩阵 $A$ 进行奇异值分解, 得到 $A=U \Sigma V^{\mathrm{T}}, V$ 和 $U$ 都是正交矩阵, 所以 $V$ 的 列向量 $v_1, v_2, \cdots, v_n$ 构成 $\mathbf{R}^n$ 空间的一组标准正交基, 表示 $\mathbf{R}^n$ 中的正交坐标系的 旋转或反射变换; $U$ 的列向量 $u_1, u_2, \cdots, u_m$ 构成 $\mathbf{R}^m$ 空间的一组标准正交基, 表示 $\mathbf{R}^m$ 中的正交坐标系的旋转或反射变换; $\Sigma$ 的对角元素 $\sigma_1, \sigma_2, \cdots, \sigma_n$ 是一组非负实 数, 表示 $\mathbf{R}^n$ 中的原始正交坐标系坐标轴的 $\sigma_1, \sigma_2, \cdots, \sigma_n$ 倍的缩放变换。

任意一个向量 $x \in \mathbf{R}^n$, 经过基于 $A=U \Sigma V^{\mathrm{T}}$ 的线性变换, 等价于经过坐标系 的旋转或反射变换 $V^{\mathrm{T}}$, 坐标轴的缩放变换 $\Sigma$, 以及坐标系的旋转或反射变换 $U$, 得 到向量 $A x \in \mathbf{R}^m$ 。图 $15.1$ 给出直观的几何解释 (见文前彩图)。原始空间的标准正交 基 (红色与黄色), 经过坐标系的旋转变换 $V^{\mathrm{T}}$ 、坐标轴的缩放变换 $\Sigma$ (黑色 $\sigma_1, \sigma_2$ ） 坐标系的旋转变换 $U$, 得到和经过线性变换 $A$ 等价的结果。

**定义**(**奇异值分解**) 矩阵的奇异值分解是指, 将一个非零的 $m \times n$ 实矩阵 $A, A \in \mathbf{R}^{m \times n}$, 表示为以下三个实矩阵乘积形式的运算 (1), 即进行矩阵的因子分解:

$$
A=U \Sigma V^{\mathrm{T}}
$$

其中 $U$ 是 $m$ 阶**标准正交矩阵** ( orthogonal matrix), $V$ 是 $n$ 阶**标准正交矩阵**, $\Sigma$ 是由降序排列的非负的对角线元素组成的 $m \times n$ 矩形对角矩阵 (rectangular diagonal matrix)

其中：

$$
\begin{aligned}
&\mathrm{UU}^{\mathrm{T}}=\mathrm{I} \\
&\mathrm{V} \mathrm{V}^{\mathrm{T}}=\mathrm{I} \\
&\Sigma=\operatorname{diag}\left(\sigma_1, \sigma_2, \cdots, \sigma_{\mathrm{p}}\right) \\
&\sigma_1 \geqslant \sigma_2 \geqslant \cdots \geqslant \sigma_{\mathrm{p}} \geqslant 0 \\
&\mathrm{p}=\min (\mathrm{m}, \mathrm{n})
\end{aligned}
$$

如下图：

<div align=center>
<img src="https://user-images.githubusercontent.com/111955215/195971524-e2cac564-c772-4e7f-b326-9f18c2e74cb4.png" width="500">
</div>

