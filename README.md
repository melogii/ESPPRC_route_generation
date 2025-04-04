# ESPPRC solver for route generation

This repository contains an implementation of the **Elementary Shortest Path Problem with Resource Constraints (ESPPRC)**, commonly used in column generation approaches for vehicle routing problems like the VRPTW (Vehicle Routing Problem with Time Windows).

This version is speacially designed for route generation in the context of the Pickup and Delivery Problem (PDP). It is meant to be used in the phase of the Branch and Price algorithm that decides which variables to include to the next instance of the problem.

## ✨ About ESPPRC

The ESPPRC consists of finding the shortest path in a graph from a source node to a destination node, while:
- Visiting each node at most once (**elementary path**),
- Respecting resource constraints (e.g., time windows, capacities),
- And minimizing the cost (e.g., distance, time, etc.).

It is a key subproblem in **Branch-and-Price** and **Column Generation** methods for solving large-scale combinatorial problems.
