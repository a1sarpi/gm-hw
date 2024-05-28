function period_progon(A, d) {
    let N = A.length;
    let p = new Array(N - 1),
        r = new Array(N - 1),
        q = new Array(N);

    let a = new Array(N), b = new Array(N),
        c = new Array(N);
    a[0] = A[0][N-1];
    c[N-1] = A[N-1][0];
    for (let i = 0; i < N; ++i) {
        b[i] = A[i][i];
        if (i != 0)
            a[i] = A[i][i-1];
        if (i != N-1)
            c[i] = A[i][i+1];
    }

    p[0] = c[0] / b[0];
    r[0] = a[0] / b[0];
    q[0] = d[0] / b[0];
    for (let i = 1; i < N; ++i) {
        p[i] = c[i] / (b[i] - a[i] * p[i-1]);
        r[i] = - a[i] * r[i-1] / (b[i] - a[i] * p[i-1]);
        q[i] = (d[i] - a[i] * q[i-1]) / (b[i] - a[i] * p[i-1]);
    }

    let s = new Array(N),
        t = new Array(N);
    s[N-1] = 1;
    t[N-1] = 0;
    for (let i = N-2; i >= 0; --i) {
        s[i] = - p[i] * s[i+1] - r[i];
        t[i] = q[i] - p[i] * t[i+1];
    }

    console.log(a, b, c, d, p, r, q, s, t);

    let x = new Array(N);
    x[N-1] = ( - c[N-1] * t[0] - a[N-1] * t[N-2] + d[N-1] ) / (c[N-1] * s[0] + a[N-1] * s[N-2] + b[N-1]);
    for (let i = 0; i < N-1; ++i) {
        x[i] = s[i] * x[N-1] + t[i];
    }
    return x;

}
