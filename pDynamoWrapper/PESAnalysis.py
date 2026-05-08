#!/usr/bin/env python3
# -*- coding: utf-8 -*-



from importlib.resources import path
import math
from multiprocessing import heap
import numpy as np

#-------------------------------------------------------------------------------
class PESAnalysis:
    '''
    '''
    def __init__(self, z, in_point, fin_point, _xlength, _ylength):
        self.z = z 
        self.xlen = _xlength
        self.ylen = _ylength
        self.start_p = in_point
        self.end_p   = fin_point
        self.structures = { "reactants": None, "products": None, "saddle_points": [] }
        self.pathx = []
        self.pathy = []

    #--------------------------------------------------------------------------
    def Analysis(self):
        '''
        '''
        self.Find_Minima()
        self.MEP()
    #--------------------------------------------------------------------------
    def Find_Minima(self):
        '''
        '''
        minima = []
        for i in range(0, self.ylen - 1):
            for j in range(0, self.xlen - 1):
                center = self.z[i, j]
                # Check 8 neighbors
                neighbors = []
                if i > 0 and j > 0:
                    neighbors = [
                        self.z[i-1, j-1], self.z[i-1, j], self.z[i-1, j+1],
                        self.z[i, j-1],                self.z[i, j+1],
                        self.z[i+1, j-1], self.z[i+1, j], self.z[i+1, j+1]
                    ]
                if i==0 and j==0:
                    neighbors = [ self.z[i, j+1], self.z[i+1, j], self.z[i+1, j+1] ]
                elif j==0:
                    neighbors = [self.z[i-1, j], self.z[i-1, j+1], self.z[i, j+1], self.z[i+1, j], self.z[i+1, j+1]]
                elif i==0:
                    neighbors = [self.z[i, j-1], self.z[i, j+1], self.z[i+1, j-1], self.z[i+1, j], self.z[i+1, j+1]]
                elif i==self.ylen-1 and j==0:
                    neighbors = [ self.z[i-1, j], self.z[i-1, j+1], self.z[i, j+1] ]
                elif i==self.ylen-1:
                    neighbors = [self.z[i-1, j-1], self.z[i-1, j], self.z[i-1, j+1], self.z[i, j-1], self.z[i, j+1]]
                elif j==self.xlen-1:
                    neighbors = [self.z[i-1, j-1], self.z[i-1, j], self.z[i, j-1], self.z[i+1, j-1], self.z[i+1, j]]                
                if center < min(neighbors):
                    minima.append((j, i, center))
        
        print(f"  Found {len(minima)} local minima")
        for x, y, e in minima[:10]:  # Print first 10
            print(f"    - Position ({x:3d}, {y:3d}): Energy = {e:10.4f} kJ/mol")
        if len(minima) > 10:
            print(f"    ... and {len(minima)-10} more")
        
        # Find reactants and products if initial/final points are provided
        if  len(minima) > 0:
            if self.start_p is None:
                self.start_p = (minima[0][0], minima[0][1])
            print(f"\n  Finding nearest minimum to initial point {self.start_p}...")
            nearest_reactant = min(minima, key=lambda m: (m[0]-self.start_p[0])**2 + (m[1]-self.start_p[1])**2)
            self.structures['reactants'] = nearest_reactant
            dist_to_reactant = np.sqrt((nearest_reactant[0]-self.start_p[0])**2 + (nearest_reactant[1]-self.start_p[1])**2)
            if dist_to_reactant > 5.0:
                print(f"    → Reactants: ({self.start_p[0]}, {self.start_p[1]}), Energy = {self.z[self.start_p[1],self.start_p[0]]:.6f} kJ, Distance = {dist_to_reactant:.2f}")
                self.structures['reactants'] = ( self.start_p[0], self.start_p[1], self.z[self.start_p[1],self.start_p[0]] )
            else:
                print(f"    → Reactants: ({nearest_reactant[0]}, {nearest_reactant[1]}), Energy = {nearest_reactant[2]:.6f} kJ, Distance = {dist_to_reactant:.2f}")
        
        
        if  len(minima) > 0:
            if self.end_p is None:
                self.end_p = (minima[-1][0], minima[-1][1])
            print(f"\n  Finding nearest minimum to final point {self.end_p}...")
            nearest_product = min(minima, key=lambda m: (m[0]-self.end_p[0])**2 + (m[1]-self.end_p[1])**2)
            self.structures['products'] = nearest_product
            dist_to_product = np.sqrt((nearest_product[0]-self.end_p[0])**2 + (nearest_product[1]-self.end_p[1])**2)
            if dist_to_product > 5.0:
                print(f"    → Products: ({self.end_p[0]}, {self.end_p[1]}), Energy = {self.z[self.end_p[1],self.end_p[0]]:.6f} kJ, Distance = {dist_to_product:.2f}")
                self.structures['products'] = ( self.end_p[0], self.end_p[1], self.z[self.end_p[1],self.end_p[0]] )
            else:
                print(f"    → Products: ({nearest_product[0]}, {nearest_product[1]}), Energy = {nearest_product[2]:.6f} kJ, Distance = {dist_to_product:.2f}")
        
    #-----------------------------------------------------------------------------------------
    def Find_Saddle_Point(self,start_xy, goal_xy):
            
        """
        Path that minimises the maximum energy (lowest barrier).
        Uses a variant of Dijkstra where the cost = max(max_so_far, energy_of_cell).
        This finds the true minimum-energy barrier path (MEP).
        
        Args:
            start_xy: (x, y) grid indices of starting point.
            goal_xy:  (x, y) grid indices of goal point.
        
        Returns:
            path: list of (x, y) indices from start to goal.
            energies: list of energy values along the path.
            barrier: maximum energy along path minus start energy (kcal/mol).
        """
        import heapq
        
        if start_xy[0] >= self.xlen or start_xy[1] >= self.ylen or goal_xy[0] >= self.xlen or goal_xy[1] >= self.ylen:
            print("Error: Start or goal coordinates are out of bounds.")
            return -1

        print("\n" + "="*70)
        print("FINDING MINIMUM BARRIER PATH (Minimax Algorithm)")
        print("="*70)
        print(f"Start: {start_xy}, Goal: {goal_xy}")
        print(f"Start energy: {self.z[start_xy[1]][start_xy[0]]:.4f} kJ/mol")
        print(f"Goal energy:  {self.z[goal_xy[1]][goal_xy[0]]:.4f} kJ/mol")
        print(f"Algorithm: Minimize maximum energy along path (true MEP)\n")
        
        rows, cols = self.ylen, self.xlen
        z = self.z
        
        print(f"Grid dimensions: {cols} x {rows}")
        
        # best_max[y][x] = lowest possible maximum energy to reach (x,y)
        best_max = np.full((rows, cols), np.inf)
        prev = np.full((rows, cols), None, dtype=object)
        
        sx, sy = start_xy[0], start_xy[1]
        best_max[sy][sx] = z[sy][sx]        # initial max = start energy
        heap = [(z[sy][sx], sx, sy)]        # (current_max, x, y)
        
        print(f"Starting Minimax search...")
        nodes_explored = 0
        nodes_added = 0
        
        while heap:
            current_max, x, y = heapq.heappop(heap)
            nodes_explored += 1
        
            # Stop if we reached the goal
            if (x, y) == goal_xy:
                print(f"\n✓ Goal reached at ({x}, {y})")
                print(f"  Minimum barrier: {current_max:.4f} kJ/mol")
                print(f"  Nodes explored: {nodes_explored}")
                break
        
            # If this entry is outdated, skip it
            if current_max > best_max[y][x]:
                continue
        
            # Periodic status output
            if nodes_explored % 100 == 0:
                print(f"  Progress: {nodes_explored} nodes explored, heap size: {len(heap)}, current position: ({x}, {y}), max barrier: {current_max:.4f}")
        
            # Explore the 8 neighbours (includes diagonals)
            for dx in (-1, 0, 1):
                for dy in (-1, 0, 1):
                    if dx == 0 and dy == 0:
                        continue
                    nx, ny = x + dx, y + dy
                
                    # Check grid boundaries
                    if 0 <= nx < cols and 0 <= ny < rows:
                    # Minimax: new barrier = max of current barrier and this cell's energy
                        new_max = max(current_max, z[ny][nx])
                    
                        # If we found a path with lower barrier to (nx, ny), update and push to heap
                        if new_max < best_max[ny][nx]:
                            best_max[ny][nx] = new_max
                            prev[ny][nx] = (x, y)
                            heapq.heappush(heap, (new_max, nx, ny))
                            nodes_added += 1
        
        print(f"  Total nodes added to heap: {nodes_added}\n")
        
        # Reconstruct path from start to goal by walking backwards
        if prev[goal_xy[1]][goal_xy[0]] is None:
            print("✗ Error: Goal not reachable from start.")
            return [], [], np.inf
        
        print(f"Reconstructing path...")
        path = []
        cur = goal_xy
        while cur is not None:
            path.append(cur)
            cur = prev[cur[1]][cur[0]]
        
        path.reverse()
        print(f"✓ Path reconstructed: {len(path)} waypoints\n")
        
        # Extract energies along the path
        energies = [z[y][x] for (x, y) in path]
        self.energies1D = energies
        
        # Calculate barrier height (kJ/mol, convert to kcal/mol for display)
        barrier = best_max[goal_xy[1]][goal_xy[0]] - z[sy][sx]
        barrier_kcal = barrier / 4.184  # Convert kJ/mol to kcal/mol
        
        # Find the barrier point (maximum energy along the path)
        barrier_idx = max(range(len(energies)), key=lambda i: energies[i])
        barrier_xy = path[barrier_idx]
        barrier_x, barrier_y = barrier_xy[0], barrier_xy[1]
        barrier_energy_max = energies[barrier_idx]
        
        print(f"Path Statistics:")
        print(f"  Path length: {len(path)} steps")
        print(f"  Start energy: {energies[0]:.6f} kJ/mol")
        print(f"  End energy: {energies[-1]:.6f} kJ/mol")
        print(f"  Max energy on path: {max(energies):.6f} kJ/mol")
        print(f"  Min energy on path: {min(energies):.6f} kJ/mol")
        print(f"  Barrier height: {barrier:.6f} kJ/mol = {barrier_kcal:.4f} kcal/mol")
        print(f"  Barrier point: ({barrier_x}, {barrier_y}) at index {barrier_idx}")
        print(f"  Barrier energy: {barrier_energy_max:.6f} kJ/mol")
        print("="*70 + "\n")
        
        if len(self.structures["saddle_points"]) > 0:
            first_e = self.z[start_xy[1]][start_xy[0]]
            if barrier_energy_max > first_e + 1.0:  # Only add if it's a significant new saddle point 
                self.structures['saddle_points'].append( (barrier_x, barrier_y, barrier_energy_max) )
            else:
                print(f"  Found barrier at ({barrier_x}, {barrier_y}) with energy {barrier_energy_max:.6f} kJ/mol, but it's not significantly higher than last saddle point ({first_e:.6f} kJ/mol). Not adding to saddle points list.")
                return -1
        else:
            self.structures['saddle_points'].append( (barrier_x, barrier_y, barrier_energy_max) )

        return path
    
    #-----------------------------------------------------------------------------------------
    def MEP(self): 
        '''
        '''
        if self.structures['reactants'] is None or self.structures['products'] is None:
            print("Error: Reactant and product structures must be identified before finding MEP.")
            return
        
        start_xy = (self.structures['reactants'][0], self.structures['reactants'][1])
        goal_xy = (self.structures['products'][0], self.structures['products'][1])
        
        #Find all saddle points along the path from reactants to products
        
        dirs = [ (1, 0), (0, 1), (1, 1) ]  # right, down, diagonal
        max_iterations = 4
        iteration = 0
        start_saddle = (start_xy[0], start_xy[1])
        if len(self.structures['saddle_points']) == 0:
            print(f"Finding first saddle point between reactants at {start_xy} and products at {goal_xy}...")
            self.Find_Saddle_Point(start_saddle, goal_xy)
        if len(self.structures['saddle_points']) > 0:
            print(f"First saddle point found at ({self.structures['saddle_points'][-1][0]}, {self.structures['saddle_points'][-1][1]}), Energy = {self.structures['saddle_points'][-1][2]:.6f} kJ/mol")
            start_saddle = ( self.structures['saddle_points'][-1][0]+1, self.structures['saddle_points'][-1][1]+1 )

            while (self.structures["saddle_points"][-1][0], self.structures["saddle_points"][-1][1]) != (goal_xy[0], goal_xy[1]):
                result = self.Find_Saddle_Point(start_saddle, goal_xy)
                if result == -1:
                    print("Failed to find next saddle point. Stopping.")
                    break
                start_saddle = ( self.structures["saddle_points"][-1][0], self.structures["saddle_points"][-1][1] )
                iteration += 1 
                if iteration >= max_iterations:
                    print("Reached maximum iterations while finding MEP. Stopping.")
                    break
        
        print(f"Total saddle points found along path: {len(self.structures['saddle_points'])}")
        print("Saddle points (x, y, energy):")
        for saddle in self.structures['saddle_points']:
            print(f"  ({saddle[0]}, {saddle[1]}, {saddle[2]})")
        # Now we have the reactants, products, and all saddle points along the path. We can construct the MEP by connecting these points.
        # For simplicity, we will just connect the points in order: reactants → saddle points → products. In a real implementation, you might want to refine the path between these points.        

        cp = [ self.structures['reactants'][0], self.structures['reactants'][1] ]
        final_points = []
        for saddle in self.structures['saddle_points']:
            final_points.append( ( saddle[0], saddle[1] ) )
        final_points.append( ( self.structures['products'][0], self.structures['products'][1] ) )

        for fin_point in final_points:
            print(final_points)
            print(f"\nFinding path from ({cp[0]}, {cp[1]}) to ({fin_point[0]}, {fin_point[1]})...")
            while not cp == fin_point:
                
                if cp[0] >= fin_point[0] and cp[1] >= fin_point[1]:
                    print( "Reached the final point: {} {}".format(cp[0], cp[1] ) )
                    print( "Energy of the final point: {}".format( self.z[ cp[1],cp[0] ] ) )
                    break

                try:
                    print( "Current Point is: {} {} ".format(cp[0], cp[1] ) )
                    print( "Energy of the Current Point: {}".format( self.z[ cp[1],cp[0] ] ) )
                except IndexError:
                    print(f"Error: Current point ({cp[0]}, {cp[1]}) is out of bounds. Stopping pathfinding.")
                    break
    
                A = math.inf
                B = math.inf
                C = math.inf            
    
                if ( cp[0] + 1 ) < self.xlen: 
                    if  (cp[0] + 1) <= fin_point[0]:
                        A = self.z[ cp[1], cp[0] ] + self.z[ cp[1], (cp[0] + 1) ]
                        print( "Increment in X:  {}".format(A) )
                else: print( "No Increment in X:  {}".format(A) )
                if ( cp[1] + 1 ) < self.ylen:
                    if  (cp[1] + 1) <= fin_point[1] : 
                        B = self.z[ cp[1], cp[0] ] + self.z[ (cp[1] + 1), cp[0] ] 
                        print( "Increment in Y:  {}".format(B) )
                else: print( "No Increment in Y: {}".format(B) )
    
                if ( cp[0] + 1 ) < self.xlen and ( cp[1] + 1 ) < self.ylen: 
                    if  (cp[0] + 1) <=fin_point[0] and (cp[1] + 1) <= fin_point[1]:
                        C = self.z[ cp[1], cp[0] ] + self.z[ (cp[1] + 1), (cp[0] + 1) ] 
                        print( "Increment in both directions:  {}".format(C) )
                else: print("No Increment in both directions: {}".format(C) )
    
                D = [ A, B, C ]
                ind = D.index(min(D))
                cp[0] = cp[0] + dirs[ind][0]
                cp[1] = cp[1] + dirs[ind][1]            
                          
                self.pathx.append(cp[0])
                self.pathy.append(cp[1])

                print(cp, fin_point)
                if cp[0] == fin_point[0] and cp[1] == fin_point[1]:
                    print( "Reached the final point: {} {}".format(cp[0], cp[1] ) )
                    print( "Energy of the final point: {}".format( self.z[ cp[1],cp[0] ] ) )
                    break
            

#-----------------------------------------------------------------------------------------
    
                
