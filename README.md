# Точная формулировка задания
1. Напишите программу, которая строит четверть поверхности цилиндра с помощью
рациональной поверхности Безье (рис. 7). Рациональная поверхность Безье строится
на базе точек, имеющих однородные координаты
$p_{00} = (1, 1, 0, 1)$, $p_{10} = (1, 1, 1, 1)$,
$p_{20} = (2, 0, 2, 2)$, $p_{01} = (−1, 1, 0, 1)$,
$p_{11} = (−1, 1, 1, 1)$, $p_{21} = (−2, 0, 2, 2)$.
Воспользуйтесь программой л.р. № 2 (Примечание: координаты
контрольных точек задаются в функции `generateControlPoints`).

2. Найдите координаты нормалей в рассчитанных точках сплайновой поверхности,
построенной в домашнем задании №1. Визуализируйте эту поверхность с
использованием заданной модели освещения.

# Генерация контрольных точек 
Происходит в файле `hw.js` в методе `generateControlPoints`, примерный код:
```javascript
this.add_coords(0, 0, Xmin, Radius, 0, 1);
this.add_coords(1, 0, Xmin, Radius, Radius, 2);
this.add_coords(2, 0, Xmin, 0, Radius, 1);

this.add_coords(0, 1, Xmax, Radius, 0, 1);
this.add_coords(1, 1, Xmax, Radius, Radius, 2);
this.add_coords(2, 1, Xmax, 0, Radius, 1);
```

Функция this.add_coords просто добавляет точку в двумерный массив,
первые два аргумента функции -- индексы i и j в этом массиве, далее x, y, z 
-- координаты точки, а последнее число -- вес точки.

Если чуть присмотреться, координаты тех точек, которые генерируются этим кодом,
несколько отличаются от тех, что даны в задании. Кажется, в задании ошибка, 
потому что точки из задания никаким образом не соотносяться с рисунком и с примерным
представлением о том, что такое цилиндр.

# Рациональный сплайн Безье

Рациональный сплайн Безье -- способ представления некоторой поверхности $ x(u, v), y(u, v), z(u, v) $.
Координаты $u, v \in [0, 1]$.

Как сказано в работе Вити, рациональный сплайн Безье представляется формулой:
$$ \mathbf{S}^{\omega} (u, v) = \dfrac{ \sum_{i=0}^n \sum_{j=0}^m B_{i, n} (u) B_{j, m} (v) \omega_{i, j} \mathbf{P}_{i, j} }
                                      { \sum_{i=0}^n \sum_{j=0}^m B_{i, n} (u) B_{j, m} (v) \omega_{i, j}}, $$
где $\mathbf{S}^\omega (u, v) = \begin{pmatrix} x(u, v) \\ y(u, v) \\ z(u, v) \end{pmatrix}$ --
координаты очередной точки (красной); $B_{k, l} (t)$ -- полиномы Бернштейна, о которых позднее;
$\mathbf{P}_{i, j} = \begin{pmatrix} x_{i, j} \\ y_{i, j} \\ z_{i, j} \end{pmatrix}$ -- 
координаты контрольных (чёрных) точек, которые как раз и генерируются с помощью `generateControlPoints`.

# Полиномы Бернштейна

Явный вид для полиномов Бернштейна:
$$ B_{i, n}(u) = C^i_n u^i (1-u)^{n-i}, $$
где $C_n^i = \dfrac{n!}{i! (n-i)!}$ -- биномиальные коэффициенты.

Подсчёт значений этого полинома реализован в функции `Bernstein` в файле `rational-bezier.js` -- 
реализован тоже как-то ускоренно, просто взял у Вити;
однако по факту используется функция `allBernstein`, которая чуть быстрее сразу считает весь вектор $B_{i, n} (u)$
для всех необходимых i.

В работе Вити предлагается использовать алгоритм де Кастельжо для подсчёта значений этих
полиномов. Этот алгоритм чуть ускоряет подсчёт биномиальных коэффициентов, но значительного прироста не
должен давать, однако где-то в коде он реализован, но не используется.

# Подсчёт точек сплайна

После того, как разобрались со всеми выражениями, входящими в формулу, можно приступать к 
отрисовке точек сплайна. Для этого необходимо разделить область параметрических координат $(u, v)$ на квадратики,
причём пусть по $u$ будет $N$ квадратиков, а по $v$ -- $M$. $N$ и $M$ вынесем в гуишку.

Тогда очередная точка сплайна с индексами $i, j$ будет иметь параметрические координаты
$u = i / N, v = j / M$. Ниже приведён кусок кода функции `calculateRationalBezierSpline`:
```javascript
const du = 1. / (N-1);
const dv = 1. / (M-1);

for (let i = 0; i < N; i++) {
    let u = i * du;
    for (let j = 0; j < M; j++) {
        let v = j * dv;

        let H = 0; // знаменатель дроби, общий для x, y, z
        let Bn = allBernstein(N_ctr-1, u);
        let Bm = allBernstein(M_ctr-1, v);
        for (let r = 0; r < N_ctr; ++r)
            for (let s = 0; s < M_ctr; ++s)
                H += Bn[r]*Bm[s]*this.pointsCtr[r][s].omega;
        let x = 0, y = 0, z = 0;
        for (let r = 0; r < N_ctr; ++r) {
            for (let s = 0; s < M_ctr; ++s) {
                x += Bn[r] * Bm[s] * this.pointsCtr[r][s].omega * this.pointsCtr[r][s].x;
                y += Bn[r] * Bm[s] * this.pointsCtr[r][s].omega * this.pointsCtr[r][s].y;
                z += Bn[r] * Bm[s] * this.pointsCtr[r][s].omega * this.pointsCtr[r][s].z;
            }
        }

        x = x/H;
        y = y/H;
        z = z/H;

        this.pointsSpline[i][j] = new Point(x, y, z);
    }
}
```

# Производная полиномов Бернштейна

Для $i \neq 0, n$:
$$ B'_{i, n}(u) = C^i_n [i u^{i-1} (1-u)^{n-i} - (n-i) u^i (1-u)^{n-i-1}] = C^i_n u^{i-1} (1-u)^{n-i-1} (i - un)$$

Заметим, что
$$ B'_{i, n} (u) = \dfrac{n}{n-i} C^i_{n-1} u^{i-1} (1-u)^{n-i-1} (i - un) = \dfrac{n (i-un)}{(n-1) u} B_{i, n-1}(u) $$
-- используем это соотношение, чтобы лишний раз не считать факториалы.

Для $i = 0$:
$$B_{0, n} (u) = (1-u)^n; \quad B'_{0, n} (u) = - n (1-u)^{n-1}$$

Для $i = n$:
$$B_{n, n} (u) = u^n; \quad B'_{0, n} (u) = n u^{n-1}.$$

Производная реализована в функции `Bernstein_prime` в файле `rational_bezier.js`.

# Нахождение касательных векторов

Для нахождения нормалей нам понадобиться найти касательные вектора
$\dfrac{\partial \mathbf{S}}{\partial u}$ и $\dfrac{\partial \mathbf{S}}{\partial v}$,
как видно, для нахождения этих производных достаточно заменить один из множителей
с обычного Бернштейна на его производную:
$$ x'_u(u, v) = \left(\dfrac{X(u, v)}{H(u, v)}\right)'_u = \dfrac{X'_u H - X H'_u}{H^2}, $$
$H$ и $X$ для данной точки мы уже считали,
$X'$ считается аналогично $X$, но вместо $B_{i, n} (u)$ подставляются $B'_{i, n} (u)$.

Код для подсчёта касательных векторов:
```javascript
let Bn_prime_u = allBernstein_prime(N_ctr-1, u);
let Bm_prime_v = allBernstein_prime(M_ctr-1, v);


let H_prime_u = 0,
    x_prime_u = 0,
    y_prime_u = 0,
    z_prime_u = 0,
    H_prime_v = 0,
    x_prime_v = 0,
    y_prime_v = 0,
    z_prime_v = 0;
for (let r = 0; r < N_ctr; ++r)
for (let s = 0; s < M_ctr; ++s) {
    H_prime_u += Bn_prime_u[r]*Bm[s]*this.pointsCtr[r][s].omega;
    H_prime_v += Bn[r]*Bm_prime_v[s]*this.pointsCtr[r][s].omega;
}
for (let r = 0; r < N_ctr; ++r) {
    for (let s = 0; s < M_ctr; ++s) {
        x_prime_u += Bn_prime_u[r] * Bm[s] * this.pointsCtr[r][s].omega * this.pointsCtr[r][s].x;
        y_prime_u += Bn_prime_u[r] * Bm[s] * this.pointsCtr[r][s].omega * this.pointsCtr[r][s].y;
        z_prime_u += Bn_prime_u[r] * Bm[s] * this.pointsCtr[r][s].omega * this.pointsCtr[r][s].z;

        x_prime_v += Bn[r] * Bm_prime_v[s] * this.pointsCtr[r][s].omega * this.pointsCtr[r][s].x;
        y_prime_v += Bn[r] * Bm_prime_v[s] * this.pointsCtr[r][s].omega * this.pointsCtr[r][s].y;
        z_prime_v += Bn[r] * Bm_prime_v[s] * this.pointsCtr[r][s].omega * this.pointsCtr[r][s].z;
    }
}

//CALCULATE TANGENT VECTORS
const x_u = (x_prime_u * H - x * H_prime_u) / (H*H);
const y_u = (y_prime_u * H - y * H_prime_u) / (H*H);
const z_u = (z_prime_u * H - z * H_prime_u) / (H*H);

const x_v = (x_prime_v * H - x * H_prime_v) / (H*H);
const y_v = (y_prime_v * H - y * H_prime_v) / (H*H);
const z_v = (z_prime_v * H - z * H_prime_v) / (H*H);

const pt_u = vec3.fromValues(x_u, y_u, z_u);
const pt_v = vec3.fromValues(x_v, y_v, z_v);
```

Код для подсчёта нормалей:
```javascript
const normal = vec3.create();
vec3.cross(normal, pt_u, pt_v);
vec3.normalize(normal, normal);

this.normalsSpline[i][j][0] = normal[0];
this.normalsSpline[i][j][1] = normal[1];
this.normalsSpline[i][j][2] = normal[2];
```
