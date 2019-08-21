---
layout: post
title: Final report for GSoC 2019 (Week 12)
comments: true
tags: [ gsoc, Sympy]
---


It’s finally the last week of the Google Summer of Code 2019. Before I start
discussing my work over the summer I would like to highlight my general
experience with the GSoC program.

GSoC gives students all over the world the opportunity to connect and
collaborate with some of the best programmers involved in open source from
around the world. I found the programme tremendusly enriching both in terms of
the depth in which I got to explore some of the areas involved in my project
and also gave me exxposure to some areas I had no previous idea about.
The role of a mentor in GSoC is the most important and I consider myself
very lucky to have got Yathartha Anirudh Joshi and Amit Kumar as my mentors.
Amit and Yathartha has been tremendously encouraging and helpful throughout the summer.
I would also like to mention the importance of the entire community involved,
just being part of the SymPy community.

### Work Completed

Here is a list of PRs which were opened during the span of GSoC:

1. [#16796 Added `_solve_modular` for handling equations a - Mod(b, c) = 0 where only b is expr](https://github.com/sympy/sympy/pull/16976)

2. [#16890 Fixing lambert in bivariate to give all real solutions](https://github.com/sympy/sympy/pull/16890)

3. [#16960 (Don't Merge)(Prototype) Adding abs while converting equation to log form to get solved by `_lambert`](https://github.com/sympy/sympy/pull/16960)

4. [#17043 Feature power_list to return all powers of a variable present in f](https://github.com/sympy/sympy/pull/17043)

5. [#17079 Defining ImageSet Union](https://github.com/sympy/sympy/pull/17079)

Here is a list of PRs merged:

1. [#16796 Added `_solve_modular` for handling equations a - Mod(b, c) = 0 where only b is expr](https://github.com/sympy/sympy/pull/16976)

2. [#16890 Fixing lambert in bivariate to give all real solutions](https://github.com/sympy/sympy/pull/16890)

Here is all the brief description about the PRs merged:

1. [#16796 Added `_solve_modular` for handling equations a - Mod(b, c) = 0 where only b is expr](https://github.com/sympy/sympy/pull/16976)

In this PR a new solver `_solve_modular` was made for solving modular equations.

### What type of equations to be considered and what domain?
```
A - Mod(B, C) = 0

    A -> This can or cannot be a function specifically(Linear, nth degree single
         Pow, a**f_x and Add and Mul) of symbol.(But currently its not a
        function of x)
    B -> This is surely a function of symbol.
    C -> It is an integer.
And domain should be a subset of S.Integers.
```
### Filtering out equations
A check is being applied named `_is_modular` which verifies that only above
mentioned type equation should return True.

### Working of `_solve_modular`
In the starting of it there is a check if domain is a subset of Integers.
```
domain.is_subset(S.Integers)
```
Only domain of integers and it subset are being considered while solving
these equations.
Now after this it separates out a modterm and the rest term on either
sides by this code.
```
modterm = list(f.atoms(Mod))[0]
rhs = -(S.One)*(f.subs(modterm, S.Zero))
if f.as_coefficients_dict()[modterm].is_negative:
    # f.as_coefficient(modterm) was returning None don't know why
    # checks if coefficient of modterm is negative in main equation.
    rhs *= -(S.One)
```
Now the equation is being inverted with the helper routine `_invert_modular`
like this.
```
n = Dummy('n', integer=True)
f_x, g_n = _invert_modular(modterm, rhs, n, symbol)
```
I am defining n in `_solve_modular` because `_invert_modular` contains
recursive calls to itself so if define the n there then it was going to have
many instances which of no use. Thats y I am defining it in `_solve_modular`.

Now after the equation is inverted now solution finding takes place.
```
if f_x is modterm and g_n is rhs:
        return unsolved_result
```
First of all if `_invert_modular` fails to invert then a ConditionSet is being
returned.
```
    if f_x is symbol:
        if domain is not S.Integers:
            return domain.intersect(g_n)
        return g_n
```
And if `_invert_modular` is fully able to invert the equation then only domain
intersection needs to takes place. `_invert_modular` inverts the equation
considering S.Integers as its default domain.

```
    if isinstance(g_n, ImageSet):
        lamda_expr = g_n.lamda.expr
        lamda_vars = g_n.lamda.variables
        base_set = g_n.base_set
        sol_set = _solveset(f_x - lamda_expr, symbol, S.Integers)
        if isinstance(sol_set, FiniteSet):
            tmp_sol = EmptySet()
            for sol in sol_set:
                tmp_sol += ImageSet(Lambda(lamda_vars, sol), base_set)
            sol_set = tmp_sol
        return domain.intersect(sol_set)
```
In this case when g_n is an ImageSet of n and f_x is not symbol so the
equation is being solved by calling `_solveset` (this will not lead to
recursion because equation to be entered is free from Mod) and then
the domain intersection takes place.

### What does `_invert_modular` do?
This function helps to convert the equation `A - Mod(B, C) = 0` to a
form (f_x, g_n).
First of all it checks the possible instances of invertible cases if not then
it returns the equation as it is.
```
a, m = modterm.args
if not isinstance(a, (Dummy, Symbol, Add, Mul, Pow)):
        return modterm, rhs
```
Now here is the check for complex arguments and returns the equation as it is
if somewhere it finds I.
```
if rhs.is_real is False or any(term.is_real is False \
            for term in list(_term_factors(a))):
        # Check for complex arguments
        return modterm, rhs
```
Now after this we check of emptyset as a solution by checking range of both
sides of equation.
As modterm can have values between [0, m - 1] and if rhs is out of this range
then emptySet is being returned.
```
if (abs(rhs) - abs(m)).is_positive or (abs(rhs) - abs(m)) is S.Zero:
        # if rhs has value greater than value of m.
        return symbol, EmptySet()
```
Now the equation haveing these types are being returned as the following
```
if a is symbol:
        return symbol, ImageSet(Lambda(n, m*n + rhs), S.Integers)

    if a.is_Add:
        # g + h = a
        g, h = a.as_independent(symbol)
        if g is not S.Zero:
            return _invert_modular(Mod(h, m), (rhs - Mod(g, m)) % m, n, symbol)

    if a.is_Mul:
        # g*h = a
        g, h = a.as_independent(symbol)
        if g is not S.One:
            return _invert_modular(Mod(h, m), (rhs*invert(g, m)) % m, n, symbol)
```
The more peculiar case is of `a.is_Pow` which is handled as following.
```
if a.is_Pow:
        # base**expo = a
        base, expo = a.args
        if expo.has(symbol) and not base.has(symbol):
            # remainder -> solution independent of n of equation.
            # m, rhs are made coprime by dividing igcd(m, rhs)
            try:
                remainder = discrete_log(m / igcd(m, rhs), rhs, a.base)
            except ValueError: # log does not exist
                return modterm, rhs
            # period -> coefficient of n in the solution and also referred as
            # the least period of expo in which it is repeats itself.
            # (a**(totient(m)) - 1) divides m. Here is link of theoram:
            # (https://en.wikipedia.org/wiki/Euler's_theorem)
            period = totient(m)
            for p in divisors(period):
                # there might a lesser period exist than totient(m).
                if pow(a.base, p, m / igcd(m, a.base)) == 1:
                    period = p
                    break
            return expo, ImageSet(Lambda(n, period*n + remainder), S.Naturals0)
        elif base.has(symbol) and not expo.has(symbol):
            remainder_list = nthroot_mod(rhs, expo, m, all_roots=True)
            if remainder_list is None:
                return symbol, EmptySet()
            g_n = EmptySet()
            for rem in remainder_list:
                g_n += ImageSet(Lambda(n, m*n + rem), S.Integers)
            return base, g_n
```
Two cases are being created based of a.is_Pow
1. x**a
2. a**x

x**a -  It is being handled by the helper function `nthroot_mod` which returns
        required solution. I am not going into very mch detail for more
        information you can read the documentation of nthroot_mod.

a**x - For this `totient` is being used in the picture whose meaning can be
       find on this [Wikipedia](https://en.wikipedia.org/wiki/Euler's_theorem)
       page. And then its divisors are being checked to find the least period
       of solutions.

2. [#16890 Fixing lambert in bivariate to give all real solutions](https://github.com/sympy/sympy/pull/16890)

This PR went through many up and downs and nearly made to the most commented PR.
And with the help of @smichr it was successfully merged. It mainly solved the
bug for not returning all solutions of lambert.

## Explaining the function `_solve_lambert` (main function to solve lambert equations)

```
Input - f, symbol, gens
OutPut - Solution of f = 0 if its lambert type expression else NotImplementedError
```
This function separates out cases as below based on the main function present in
the main equation.
```
For the first ones:
1a1) B**B = R != 0 (when 0, there is only a solution if the base is 0,
                   but if it is, the exp is 0 and 0**0=1
                   comes back as B*log(B) = log(R)
1a2) B*(a + b*log(B))**p = R or with monomial expanded or with whole
                            thing expanded comes back unchanged
   log(B) + p*log(a + b*log(B)) = log(R)
   lhs is Mul:
       expand log of both sides to give:
       log(B) + log(log(B)) = log(log(R))
1b) d*log(a*B + b) + c*B = R
   lhs is Add:
       isolate c*B and expand log of both sides:
       log(c) + log(B) = log(R - d*log(a*B + b))
```
If the equation are of type 1a1, 1a2 and 1b then the mainlog of the equation is
taken into concern as the deciding factor lies in the main logarithmic term of equation.
```
For the next two,
   collect on main exp
   2a) (b*B + c)*exp(d*B + g) = R
       lhs is mul:
           log to give
           log(b*B + c) + d*B = log(R) - g
   2b) -b*B + g*exp(d*B + h) = R
       lhs is add:
           add b*B
           log and rearrange
           log(R + b*B) - d*B = log(g) + h
```
If the equation are of type 2a and 2b then the mainexp of the equation is
taken into concern as the deciding factor lies in the main exponential term of equation.
```
3) d*p**(a*B + b) + c*B = R
   collect on main pow
   log(R - c*B) - a*B*log(p) = log(d) + b*log(p)
```
If the equation are of type 3 then the mainpow of the equation is
taken into concern as the deciding factor lies in the main power term of equation.

Eventually from all of the three cases the equation is meant to be converted to this form:-
```
f(x, a..f) = a*log(b*X + c) + d*X - f = 0 which has the
solution,  X = -c/b + (a/d)*W(d/(a*b)*exp(c*d/a/b)*exp(f/a)).
```
And the solution calculation process is done by `_lambert` function.

Everything seems flawless?? You might be thinking no modification is required. Lets
see what loopholes are there in it.

## What does PR [#16890](https://github.com/sympy/sympy/pull/16890) do?

There are basically two flaws present with the this approach.
1. Not considering all branches of equation while taking log both sides.
2. Calculation of roots should consider all roots in case having rational power.

### 1. Not considering all branches of equation while taking log both sides.

Let us consider this equation to be solved by `_solve_lambert` function.
```
-1/x**2 + exp(x/2)/2 = 0
```
So what the old `_solve_lambert` do is to convert this equation to following.
```
2*log(x) + x/2 = 0
```
and calculates its roots from `_lambert`.
But it missed this branch of equation while taking log on main equation.
```
2*log(-x) + x/2 = 0
```
Yeah you can reproduce the original equation from this equation.So basically the problem
was that it missed the branches of equation while taking log. And when does the
main equation have more than one branch?? The terms having even powers of variable x
leads to two different branches of equation.

So how it is solved?
What I has done is that before actually gets into solving I preprocess the main equation
and if it has more than one branches of equation while converting taking log then I consider
all the equations generated from them.(with the help of `_solve_even_degree_expr`)

How I preprocess the equation?
So what I do is I replace all the even powers of x present with even powers of t(dummy variable).
```
Code for targeted replacement
lhs = lhs.replace(
            lambda i:  # find symbol**even
                i.is_Pow and i.base == symbol and i.exp.is_even,
            lambda i:  # replace t**even
                t**i.exp)
Example:-
Main equation -> -1/x**2 + exp(x/2)/2 = 0
After replacement -> -1/t**2 + exp(x/2)/2 = 0
```
Now I take logarithms on both sides and simplify it.
```
After simplifying -> 2*log(t) + x/2 = 0
```
Now I call function `_solve_even_degree_expr` to replace the t with +/-x to generate two equations.
```
Replacing t with +/-x
1. 2*log(x) + x/2 = 0
2. 2*log(-x) + x/2 = 0
```
And consider the solutions of both of the equations to return all lambert real solutions
of `-1/x**2 + exp(x/2)/2 = 0`.

Hope you could understand the logic behind this work.

### 2. Calculation of roots should consider all roots in case having rational power.

This flaw is in the calculation of roots in function `_lambert`.
Earlier the function_lambert has the working like :-

1. Find all the values of a, b, c, d, e in the required loagrithmic equation
2. Then it defines a solution of the form
```
-c/b + (a/d)*l where l = LambertW(d/(a*b)*exp(c*d/a/b)*exp(-f/a), k)
```
and then it included that solution.
I agree everything seems flawless here. but try to see the step where we are defining l.

Let us suppose a hypothetical algorithm just like algorithm used in `_lambert`
in which equation to be solved is
```
x**3 - 1 = 0
```
and in which we define solution of the form
```
x = exp(I*2*pi/n) where n is the power of x in equation
```
so the algorithm will give solution
```
x = exp(I*2*pi/3) # but expected was [1, exp(I*2*pi/3), exp(-I*2*pi/3)]
```
which can be found by finding all solutions of
```
x**n - exp(2*I*pi) = 0
```
by a different correct algorithm. Thats y it was wrong.
The above algorithm would have given correct values for `x - 1 = 0`.

And the question in your mind may arise that why only exp() because the
possiblity of having more than one roots is in exp(), because if the algorithm
would have been like `x = a`, where a is some real constant then there is not
any possiblity of further roots rather than solution like `x = a**(1/n)`.
And its been done in code like this:
```
code
num, den = ((c*d-b*f)/a/b).as_numer_denom()
p, den = den.as_coeff_Mul()
e = exp(num/den)
t = Dummy('t')
args = [d/(a*b)*t for t in roots(t**p - e, t).keys()]
```

### Work under development

- [#17079 Defining ImageSet Union](https://github.com/sympy/sympy/pull/17079)

This PR tends to define a unifying algorithm for linear relations.

### Future Work
Here is a list that comprises of all the ideas (which were a part of my GSoC
Proposal and/or thought over during the SoC) which can extend my GSoC project.

1. Integrating helper solvers within solveset: linsolve, solve_decomposition, nonlinsolve

2. Handle nested trigonometric equations.
