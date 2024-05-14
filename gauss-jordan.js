function zeros(M, n, m) {
    for(let i = 0; i < n; ++i)
        for (let j = 0; j < m; ++j)
            M[i][j] = 0.;
}

function diagonalize(M) {
    let m = M.length;
    let n = M[0].length;
    for(let k=0; k<Math.min(m,n); ++k) {
        // Find the k-th pivot
        i_max = findPivot(M, k);
        if (M[i_max, k] == 0)
            throw "matrix is singular";
        swap_rows(M, k, i_max);
        // Do for all rows below pivot
        for(let i=k+1; i<m; ++i) {
            // Do for all remaining elements in current row:
            let c = M[i][k] / M[k][k];
            for(let j=k+1; j<n; ++j) {
                M[i][j] = M[i][j] - M[k][j] * c;
            }
            // Fill lower triangular matrix with zeros
            M[i][k] = 0;
        }
    }
}

function findPivot(M, k) {
    let i_max = k;
    for(let i=k+1; i<M.length; ++i)
        if (Math.abs(M[i][k]) > Math.abs(M[i_max][k]))
            i_max = i;
    return i_max;
}

function swap_rows(M, i_max, k) {
    if (i_max != k) {
        let temp = M[i_max];
        M[i_max] = M[k];
        M[k] = temp;
    }
}

function makeM(A, b) {
    let M = new Array(A.length);
    for (let i = 0; i < A.length; ++i)
        M[i] = new Array(A[i].length + 1);
    for (let i = 0; i < A.length; ++i) {
        for (let j = 0; j < A[i].length; ++j)
            M[i][j] = A[i][j];
        M[i][A[i].length] = b[i];
    }
    return M;
}

function substitute(M) {
  let m = M.length;
  for(let i=m-1; i>=0; --i) {
    let x = M[i][m] / M[i][i];
    for(let j=i-1; j>=0; --j) {
      M[j][m] -= x * M[j][i];
      M[j][i] = 0;
    }
    M[i][m] = x;
    M[i][i] = 1;
  }
}

function extractX(M) {
  let x = [];
  let m = M.length;
  let n = M[0].length;
  for(let i=0; i<m; ++i){
    x.push(M[i][n-1]);
  }
  return x;
}

function solve(A, b) {
    //print(A, "A");
    let M = makeM(A,b);
    //print(A, "M");
    diagonalize(M);
    //print(A, "diag");
    substitute(M);
    //print(A, "subst");
    let x = extractX(M);
    //print(x, "x");
    return x;
}
