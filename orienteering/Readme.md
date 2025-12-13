## The Orienteering Problem
The Mixed Integer Programming (MIP) model for the Orienteering Problem (Group 1, Case 1) is implemented in Julia and using Gurobi. Three test instances are provided in the *instances/* directory.

## Problem Description

- A set of locations $L$ with scores $s_i$ for each location $i \in L$.
- A travel time matrix where $d_{ij}$ indicates the time required to travel from location $i$ to location $j$.
- A total available time $T_{\max}$ for the route.
- A route starts and ends at depot location $1$.

### Objective
Find a feasible route that respects the time limit and maximizes the total score of the visited locations.

## Usage
Run from the command line: `julia orienteering.jl <instance-file>`

## MIP model
### Decision variables
- Binary arc variables:
$x_{ij}$ for all $i,j \in L$ with $x_{ij}=1$ if arc $(i,j)$ is used and $x_{ij}=0$ otherwise.

- Binary location variables:
$y_i$ for all $i \in L$ with $y_i=1$ if location $i$ is visited, and $y_i=0$ otherwise.

### Objective function
Maximize the collected scores:
$$\max \sum_{i \in L} s_i \cdot y_i$$

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
$$\sum_{i \in S, j \in S'} x_{ij} \geq y_k \quad \forall k \in S',\quad S'=L-S,\quad S' \neq \emptyset,\quad S \subset L \quad\text{with}\quad 1 \in S$$
$$\sum_{i \in S, j \in S'} x_{ji} \geq y_k \quad \forall k \in S',\quad S'=L-S,\quad S' \neq \emptyset,\quad S \subset L,\quad\text{with}\quad 1 \in S$$

Please note that violated subtour elimination constraints are inserted step by step using lazy cuts (isolated components are identified).

### Branching strategy
Assign higher priorities to the the $y$-variables related to the inclusion or exclusion of locations.
