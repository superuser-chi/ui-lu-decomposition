<template>
  <q-page class="row justify-center wrap">
    <div class="q-pa-md">
      <div class="text-h3 text-bold tex-italic">Introduction</div>
      <p class="text">
        In numerical analysis and linear algebra, LU decomposition (where ”LU” stands for ”lower–upper”)
        factors a matrix as the product of a lower triangular matrix and an upper triangular matrix.
        Once the matrices have been cast into LU form, they become easy to invert. The LU decomposition
        is used for solving systems of linear equations. In this section, we consider three
        algorithms that are used to decompose matrices into LU form.
      </p>
      <div class="text-h4 text-bold tex-italic">Doolittle Factorisation</div>
      <p class="text">
        If the lower triangular matrix L has 1’s on its main diagonal, then the LU factorisation is
        called a Doolittle factorization. The structure of a Doolittle factorisation of a matrix A of size
        n × n is shown below:
        $$
        \begin{equation*}
        \begin{pmatrix}
        a_{1,1} & a_{1,2} & a_{1,3} & \cdots & a_{1,n}\\
        a_{2,1} & a_{2,2} & a_{2,3} & \cdots & a_{2,n}\\
        a_{3,1} & a_{3,2} & a_{3,3} & \cdots & a_{3,n}\\
        \cdots & \cdots & \cdots & \cdots & \cdots\\
        a_{n,1} & a_{n,2} & a_{n,3} & \cdots & a_{n,n}\\
        \end{pmatrix}
        =
        \begin{pmatrix}
        1 & 0 & 0 & \cdots & 0 \\
        l_{2,1} & 1 & 0 & \cdots & 0 \\
        l_{3,1} & l_{3,2} & 1 & \cdots & 0 \\
        \cdots & \cdots & \cdots & \ddots & \vdots \\
        l_{1,1} & l_{1,2} & l_{1,3} & \cdots & 1 \\
        \end{pmatrix}
        \begin{pmatrix}
        u_{1,1} & u_{1,2} & u_{1,3} & \cdots & u_{1,n}\\
        0 & u_{2,2} & u_{2,3} & \cdots & u_{2,n}\\
        0 & 0 & u_{3,3} & \cdots & u_{3,n}\\
        \cdots & \cdots & \cdots & \ddots & \vdots\\
        0 & 0 & 0 & \cdots & u_{n,n}\\
        \end{pmatrix}
        \end{equation*}
        $$
      </p>
      <div class="text-h4 text-bold tex-italic">Crout Factorisation</div>
      <p class="text">
        If the upper triangular matrix U has 1’s on its main diagonal, then the LU factorisation is
        called a Crout factorisation. The Crout factorisation structure is shown below:
        $$
        \begin{equation*}
        \begin{pmatrix}
        a_{1,1} & a_{1,2} & a_{1,3} & \cdots & a_{1,n}\\
        a_{2,1} & a_{2,2} & a_{2,3} & \cdots & a_{2,n}\\
        a_{3,1} & a_{3,2} & a_{3,3} & \cdots & a_{3,n}\\
        \cdots & \cdots & \cdots & \cdots & \cdots\\
        a_{n,1} & a_{n,2} & a_{n,3} & \cdots & a_{n,n}\\
        \end{pmatrix}
        =
        \begin{pmatrix}
        l_{1,1} & 0 & 0 & \cdots & 0 \\
        l_{2,1} & l_{2,2} & 0 & \cdots & 0 \\
        l_{3,1} & l_{3,2} & l_{3,3} & \cdots & 0 \\
        \cdots & \cdots & \cdots & \ddots & \vdots \\
        l_{1,1} & l_{1,2} & l_{1,3} & \cdots & l_{n,n} \\
        \end{pmatrix}
        \begin{pmatrix}
        1 & u_{1,2} & u_{1,3} & \cdots & u_{1,n}\\
        0 & 1 & u_{2,3} & \cdots & u_{2,n}\\
        0 & 0 & 1 & \cdots & u_{3,n}\\
        \cdots & \cdots & \cdots & \ddots & \vdots\\
        0 & 0 & 0 & \cdots & 1\\
        \end{pmatrix}
        \end{equation*}
        $$
      </p>
      <div class="text-h4 text-bold tex-italic">Naive Gaussian Elimination</div>
      <p class="text">
        Consider the problem of solving the matrix equation $AX = b$
        $$
        \begin{equation*}
        \begin{bmatrix}
        2 & 1 & -1 & 2\\
        4 & 5 & -3 & 6\\
        -2 & 5 & -2 & 6\\
        4 & 11 & -4 & 8\\
        \end{bmatrix}
        \begin{bmatrix}
        x_1 \\
        x_2 \\
        x_3 \\
        x_4 \\
        \end{bmatrix}
        =
        \begin{bmatrix}
        -3 \\
        -7 \\
        -10 \\
        -2 \\
        \end{bmatrix}
        \end{equation*}
        $$
        A series of elementary row operations (known as forward elimination) can reduce the
        coefficient matrix to upper triangular form
        $$
        \begin{equation*}
        \begin{bmatrix}
        2 & 1 & -1 & 2\\
        0 & 5 & -3 & 6\\
        0 & 0 & -2 & 6\\
        0& 0 & 0 & 8\\
        \end{bmatrix}
        \begin{bmatrix}
        x_1 \\
        x_2 \\
        x_3 \\
        x_4 \\
        \end{bmatrix}
        =
        \begin{bmatrix}
        -3 \\
        -7 \\
        -10 \\
        -2 \\
        \end{bmatrix}
        \end{equation*}
        $$
        The forward elimination phase can be interpreted as multiplying the matrix equation,
        $$Ax = b$$ by a matrix $M$
        $$M(Ax) = Mb$$
        where M is selected in such a way that MA is upper triangular. That is
        $$
        \begin{equation*}
        MA
        =
        \begin{bmatrix}
        2 & 1 & -1 & 2\\
        0 & 5 & -3 & 6\\
        0 & 0 & -2 & 6\\
        0& 0 & 0 & 8\\
        \end{bmatrix}
        =
        U
        \end{equation*}
        $$
        The first step of the forward elimination procedure of the naive Gaussian elimination process gives
        $$
        \begin{equation*}
        \begin{bmatrix}
        2 & 1 & -1 & 2\\
        0 & 3 & -1 & 2\\
        0 & 6 & -3 & 8\\
        0 & 9 & 2 & 4\\
        \end{bmatrix}
        \begin{bmatrix}
        x_1 \\
        x_2 \\
        x_3 \\
        x_4 \\
        \end{bmatrix}
        =
        \begin{bmatrix}
        -3 \\
        -1 \\
        -13 \\
        -4 \\
        \end{bmatrix}
        \end{equation*}
        $$
        This result can be obtained by multiplying by a lower triangular matrix M1:
        $$ M_{1}Ax = M_{1}b$$
        where
        $$
        \begin{equation*}
        M_1
        =
        \begin{bmatrix}
        1 & 0 & 0 & 0\\
        -2 & 1 & 0 & 0\\
        1 & 0 & 1 & 0\\
        -2 & 0 & 0 & 1\\
        \end{bmatrix}
        \end{equation*}
        $$
        Note that row operations used in the first step are
        $$ R_2 → R_2 − 2R_1 $$
        $$ R_3 → R_3 + R_1 $$
        $$ R_4 → R_4 − 2R_1 $$
        The second forward elimination step gives
        $$
        \begin{equation*}
        \begin{bmatrix}
        2 & 1 & −1 & 2 \\
        0 & 3 & −1 & 2 \\
        0 & 0 & −1 & 4 \\
        0 & 0 & 1 −& 2 \\
        \end{bmatrix}
        \begin{bmatrix}
        x_1\\
        x_2\\
        x_3\\
        x_4\\
        \end{bmatrix}
        =
        \begin{bmatrix}
        −3 \\
        −1\\
        −11\\
        7\\
        \end{bmatrix}
        \end{equation*}
        $$
        This is equivalent to
        $$ M_2M_1Ax = M_2M_1b $$
        Note that the row operations used in this step are
        $$
        R_3 → R_3 − 2R_2 \\
        R_4 → R_4 − 3R_2 \\
        $$
        Thus,
        $$
        \begin{equation*}
        M_2 =
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & −2 & 1 & 0 \\
        0 & −3 & 0 & 1 \\
        \end{bmatrix}
        \end{equation*}
        $$
        The third step gives the upper triangular form
        $$
        \begin{equation*}
        \begin{bmatrix}
        2 & 1 & −1 & 2 \\
        0 & 3 & −1 & 2 \\
        0 & 0 & −1 & 4 \\
        0 & 0 & 0 & 2 \\
        \end{bmatrix}
        \begin{bmatrix}
        x_1 \\
        x_2 \\
        x_3 \\
        x_4 \\
        \end{bmatrix}
        =
        \begin{bmatrix}
        −3 \\
        −1 \\
        −11 \\
        −4 \\
        \end{bmatrix}
        \end{equation*}
        $$
        which is equivalent to
        $$ M_3M_2M_1Ax = M_3M_2M_1b .$$
        Since the row operation used is
        $$ R_4 → R_4 + R_3 $$
        It follows that
        $$
        \begin{equation*}
        M_3 =
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & 1 & 0 \\
        0 & 0 & 1 & 1 \\
        \end{bmatrix}
        \end{equation*}
        $$
        Thus, the matrix M can be obtained as a product of the three multiplier matrices:
        $$ M = M_3M_2M_1 $$
        The forward elimination process can thus be interpreted as follows:
        $ MA = U $
        $$
        \begin{align*}
        A &= M^{−1}U \\
        &= M_1^{−1}M_2^{−1}M_3^{−1}U \\
        &= LU \\
        \end{align*}
        $$
        where $ L = M_1^{−1}M_2^{−1}M_3^{−1}$. The inverses $M_1^{−1},M_2^{−1},M_3^{−1}$ can be obtained by simply changing
        the signs of the multipliers. Thus, we have
        $$
        \begin{equation*}
        L =
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        2 & 1 & 0 & 0 \\
        −1 & 0 & 1 & 0 \\
        2 & 0 & 0 & 1 \\
        \end{bmatrix}
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 2 & 1 & 0 \\
        0 & 3 & 0 & 1 \\
        \end{bmatrix}
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & 1 & 0 \\
        0 & 0 & −1 & 1 \\
        \end{bmatrix}
        =
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        2 & 1 & 0 & 0 \\
        −1 & 2 & 1 & 0 \\
        2 & 3 & −1 & 1 \\
        \end{bmatrix}
        \end{equation*}
        $$
        Given any matrix A, the naive Gaussian elimination method can be used to decompose
        the matrix to LU form as follows:
      </p>
      <div class="text-h4 text-bold tex-italic mt-5">step 1</div>
      $$
      \begin{array}{cc}
      R_2 → R_2 − 2R_1 \\
      R_3 → R_3 + R_1 \\
      R_4 → R_4 − 2R_1 \\
      \end{array}
      \begin{bmatrix}
      2 & 1 & −1 & 2 \\
      0 & 3 & −1 & 2 \\
      0 & 6 & −3 & 8 \\
      0 & 9 & −2 & 4 \\
      \end{bmatrix}
      $$
      $$
      M_1 =
      \begin{bmatrix}
      1 & 0 & 0 & 0 \\
      −2 & 1 & 0 & 0 \\
      1 & 0 & 1 & 0 \\
      −2 & 0 & 0 & 1 \\
      \end{bmatrix}
      ,M_1^{−1} =
      \begin{bmatrix}
      1 & 0 & 0 & 0 \\
      2 & 1 & 0 & 0 \\
      −1 & 0 & 1 & 0 \\
      2 & 0 & 0 & 1 \\
      \end{bmatrix}
      $$
      <div class="text-h4 text-bold tex-italic mt-5">step 2</div>
      $$
      \begin{array}{cc}
      R_3 → R_3 − 2R_2 \\
      R_4 → R_4 − 3R_2 \\
      \end{array}
      \begin{bmatrix}
      2 & 1 & −1 & 2 \\
      0 & 3 & −1 & 2 \\
      0 & 0 & −1 & 4 \\
      0 & 0 & 1 & −2 \\
      \end{bmatrix}
      $$
      $$
      M2 =
      \begin{bmatrix}
      1 & 0 & 0 & 0 \\
      0 & 1 & 0 & 0 \\
      0 & −2 & 1 & 0 \\
      0 & −3 & 0 & 1 \\
      \end{bmatrix}
      , M_2^{−1} =
      \begin{bmatrix}
      1 & 0 & 0 & 0 \\
      0 & 1 & 0 & 0 \\
      0 & 2 & 1 & 0 \\
      0 & 3 & 0 & 1 \\
      \end{bmatrix}
      $$
      <div class="text-h4 text-bold tex-italic mt-5">step 3</div>
      <p class="text">
        $$
        \begin{array}{cc}
        R_4 → R_4 + R_3 \\
        \end{array}
        \begin{bmatrix}
        2 & 1 & −1 & 2 \\
        0 & 3 & −1 & 2 \\
        0 & 0 & −1 & 4 \\
        0 & 0 & 0 & 2 \\
        \end{bmatrix}
        $$
        $$
        M_3 =
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & 1 & 0 \\
        0 & 0 & 1 & 1 \\
        \end{bmatrix}
        , M_3^{−1} =
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & 1 & 0 \\
        0 & 0 & −1 & 1 \\
        \end{bmatrix}
        $$
        $$
        U =
        \begin{bmatrix}
        2 & 1 & −1 & 2 \\
        0 & 3 & −1 & 2 \\
        0 & 0 & −1 & 4 \\
        0 & 0 & 0 & 2 \\
        \end{bmatrix}
        $$
        $$
        \begin{align}
        L &= M_1^{-1}M_2^{−1}M_3^{−1} \\
        &=
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        2 & 1 & 0 & 0 \\
        −1 & 0 & 1 & 0 \\
        2 & 0 & 0 & 1 \\
        \end{bmatrix}
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 2 & 1 & 0 \\
        0 & 3 & 0 & 1 \\
        \end{bmatrix}
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        0 & 1 & 0 & 0 \\
        0 & 0 & 1 & 0 \\
        0 & 0 & −1 & 1 \\
        \end{bmatrix} \\
        &=
        \begin{bmatrix}
        1 & 0 & 0 & 0 \\
        2 & 1 & 0 & 0 \\
        −1 & 2 & 1 & 0 \\
        2 & 3 & −1 & 1 \\
        \end{bmatrix}
        \end{align}
        $$
      </p>
      <div class="text-h2 text-bold tex-italic">Solving linear equations using LU factorisation</div>
      <p class="text mb-4">
        Assume that $A$ has a $LU$ factorization. The solution $X$ to the linear system $AX = b $, as
        follows:
        $$
        \begin{align}
        AX &= b \\
        LUX &= b \\
        LZ &= b, where Z = UX;
        \end{align}
        $$
        To solve the matrix equation for $X$, first we solve $LZ = b$ for $Z$ using the forward-substitution
        method, then solve $UX = Z$ for $X$ by backward-substitution method.
      </p>
      <q-pagination v-model="current" :max="5" :direction-links="true"></q-pagination>
    </div>
  </q-page>
</template>
<script>
export default {
  name: "Notes",
  data() {
    return {
      current: 1
    };
  },
  methods: {
    onLoad(index, done) {
      this.update_element("index");
      setTimeout(() => {
        if (this.items) {
          this.items.push({}, {}, {}, {}, {}, {}, {});
          done();
        }
      }, 2000);
    },
    update_element(element) {
      this.$nextTick(() => {
        window.MathJax.Hub.Queue([
          "Typeset",
          window.MathJax.Hub,
          this.$refs[`${element}`]
        ]);
      });
    }
  }
};
</script>
<style scoped>
section.preview {
  border: 1px solid #e0e0e0;
  padding: 15px;
  overflow-x: scroll;
}
section.preview ul {
  list-style: none;
  padding-left: 0;
  margin-left: 0;
}
section.preview ul strong {
  font-weight: 900;
  font-family: "Times New Roman", Times, serif;
  font-size: 1.25em !important;
}
</style>