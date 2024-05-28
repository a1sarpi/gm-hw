"use strict";


/* def Bernstein(i, n, u):
    temp = [0.0 for j in range(0, n+1)]
    temp[n-i] = 1.0
    u1 = 1.0-u
    for k in range(1, n+1):
        for j in range(n, k-1, -1):
            temp[j] = u1*temp[j]+u*temp[j-1]
    return temp[n]

def allBernstein(n, u):
    B = [0.0 for j in range(0, n+1)]
    B[0] = 1.0
    u1 = 1.0 - u
    for j in range(1, n+1):
        saved = 0.0
        for k in range(0, j):
            temp = B[k]
            B[k] = saved + u1*temp
            saved = u*temp
        B[j] = saved
    return B

def deCasteljau1(P,n,u):
    Q = copy.deepcopy(P)
    for k in range(1, n+1):
        for i in  range(0, n-k+1):
            Q[i] = (1.0-u)*Q[i]+u*Q[i+1]
    return Q[0]

def deCasteljau2(P, n, m, u0, v0):
    Q = []
    if(n <= m):
        for j in range(0, m+1):
            Q.append(deCasteljau1(P[j], n, u0))
        return deCasteljau1(Q, m, v0)
    else:
        for i in range(0, n+1):
            Q.append(deCasteljau1(P[:, i], m, v0))
        return deCasteljau1(Q, n, u0)
*/

function Bernstein(i, n, u) {
    let tmp = new Array(n+1);
    for (let i = 0; i < n+1; ++i)
        tmp[i] = 0.0;
    tmp[n-i] = 1.0;
    let u1 = 1.0 - u;
    for (let k = 1; k < n+1; ++k) {
        for (let j = n; j > k-1; --j) {
            tmp[j] = u1 * tmp[j] + u * tmp[j-1];
        }
    }
    return tmp[n];
}

function allBernstein(n, u) {
    let B = new Array(n+1);
    for (let i = 0; i < n+1; ++i)
        B[i] = 0.0;
    B[0] = 1.0;
    let u1 = 1.0 - u;
    for (let j = 1; j < n+1; ++j) {
        let saved = 0.0;
        for (let k = 0; k < j; ++k) {
            let tmp = B[k];
            B[k] = saved + u1 * tmp;
            saved = u*tmp;
        }
        B[j] = saved;
    }
    return B;
}

function Bernstein_prime(i, n, u) {
    if (i == 0)
        return - n * Math.pow(1 - u, n-1);
    else if (i == n) 
        return n * Math.pow(u, n-1);
    
    if (Math.abs(u) > 1e-6 && Math.abs(u) > 1e-6)
        return n * (i - u*n) * Bernstein(i, n-1, u) / ( (n-1) * u )
    else
        return 0.;
}

function allBernstein_prime(n, u) {
    let B = new Array(n+1);
    for (let i = 0; i < n+1; ++i)
        B[i] = Bernstein_prime(i, n, u);
    return B;
}

function deCasteljau1(P, n, u) {
    let Q = P.slice();
    for (let k = 1; k < n+1; ++k)
        for (let i = 0; i < n-k+1; ++i) {
            Q[i].x = (1.0-u)*Q[i].x + u*Q[i+1].x;
            Q[i].y = (1.0-u)*Q[i].y + u*Q[i+1].y;
            Q[i].z = (1.0-u)*Q[i].z + u*Q[i+1].z;
        }
    return Q[0]
}

function deCasteljau2(P, n, m, u0, v0) {
    let Q = []
    if(m <= n) {
        for (let i = 0; i < n+1; ++i)
            Q.push(deCasteljau1(P[i], m, u0))
        return deCasteljau1(Q, n, v0)
    } else {
        for (let j = 0; j < m+1; ++j)
            Q.push(deCasteljau1(P.slice(0, j), n, v0))
        return deCasteljau1(Q, m, u0)
    }
}
