## The Orienteering Problem
The Mixed Integer Programming (MIP) model for the Orienteering Problem (Group 1, Case 1) is implemented in Julia. Three sample test instances are provided in the *instances/* directory.

Run the program from the command line: `julia orienteering.jl <instance-file>`

The MIP model used is defined as follows:

$x_1$

**The Cauchy-Schwarz Inequality**\
$$x_{ij} \in \{ 0,1 \} \forall i,j \in L$$ with $x_{ij}=1$ if arc $(i,j)$ is used and $x_{ij}=0$ otherwise.

### Decision variables
- Binary arc variables:
$$
x_{ij} =
\begin{cases}
1 & \text{if arc } (i,j) \text{ is used}, \\
0 & \text{otherwise}.
\end{cases}
\quad \forall i,j \in L
$$

- Binary location variables:
\[
y_i =
\begin{cases}
1 & \text{if location } i \text{ is visited}, \\
0 & \text{otherwise}.
\end{cases}
\quad \forall i \in L
\]

### Objective function
Maximize the collected scores:
$$\max \sum_{i \in L} s_i \cdot y_i$$

### Constraints
1. **Linking arc and location variables (outgoing degree):**
$$\sum_{\substack{j \in L \\ j \neq i}} x_{ij} = y_i \quad \forall i \in L$$

2. **Linking arc and location variables (incoming degree):**
\[
\sum_{\substack{j \in L \\ j \neq i}} x_{ji} = y_i \quad \forall i \in L
\]

3. **Depot must always be visited:**
\[
y_1 = 1
\]

4. **Time limit constraint:**
\[
\sum_{i \in L} \sum_{j \in L} d_{ij} \cdot x_{ij} \;\leq\; T_{\max}
\]
