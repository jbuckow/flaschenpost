## The Orienteering Problem
The Mixed Integer Programming (MIP) model for the Orienteering Problem (Group 1, Case 1) is implemented in Julia. Three sample test instances are provided in the *instances/* directory.

Run the program from the command line: `julia orienteering.jl <instance-file>`

The MIP model used is defined as follows:

### Decision variables
- Binary arc variables:
$x_{ij}$ for all $i,j \in L$ with $x_{ij}=1$ if arc $(i,j)$ is used and $x_{ij}=0$ otherwise.

- Binary location variables:
$y_i$ for all $i \in L$ with $y_i=1$ if location $i$ is visited, and $y_i=0$ otherwise.

### Objective function
Maximize the collected scores:
$$\max \sum\limits_{i \in L} {s_i \cdot y_i}$$

### Constraints
1. **Linking arc and location variables (outgoing degree):**
$$\sum_{{j \in L \\ j \neq i}} x_{ij} = y_i \quad \forall i \in L$$

2. **Linking arc and location variables (incoming degree):**
$$\sum_{\substack{j \in L \\ j \neq i}} x_{ji} = y_i \quad \forall i \in L$$

3. **Depot must always be visited:**
$$y_1 = 1$$

4. **Time limit must be respected:**
$$\sum_{i \in L} \sum_{j \in L} d_{ij} \cdot x_{ij} \leq T_{\max}$$

5. **Subtour elimination constraints:**
$$xxx$$
