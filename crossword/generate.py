import sys
import copy
import operator

from crossword import *


class CrosswordCreator():

    def __init__(self, crossword):
        """
        Create new CSP crossword generate.
        """
        self.crossword = crossword
        self.domains = {
            var: self.crossword.words.copy()
            for var in self.crossword.variables
        }

    def letter_grid(self, assignment):
        """
        Return 2D array representing a given assignment.
        """
        letters = [
            [None for _ in range(self.crossword.width)]
            for _ in range(self.crossword.height)
        ]
        for variable, word in assignment.items():
            direction = variable.direction
            for k in range(len(word)):
                i = variable.i + (k if direction == Variable.DOWN else 0)
                j = variable.j + (k if direction == Variable.ACROSS else 0)
                letters[i][j] = word[k]
        return letters

    def print(self, assignment):
        """
        Print crossword assignment to the terminal.
        """
        letters = self.letter_grid(assignment)
        for i in range(self.crossword.height):
            for j in range(self.crossword.width):
                if self.crossword.structure[i][j]:
                    print(letters[i][j] or " ", end="")
                else:
                    print("█", end="")
            print()

    def save(self, assignment, filename):
        """
        Save crossword assignment to an image file.
        """
        from PIL import Image, ImageDraw, ImageFont
        cell_size = 100
        cell_border = 2
        interior_size = cell_size - 2 * cell_border
        letters = self.letter_grid(assignment)

        # Create a blank canvas
        img = Image.new(
            "RGBA",
            (self.crossword.width * cell_size,
             self.crossword.height * cell_size),
            "black"
        )
        font = ImageFont.truetype("assets/fonts/OpenSans-Regular.ttf", 80)
        draw = ImageDraw.Draw(img)

        for i in range(self.crossword.height):
            for j in range(self.crossword.width):

                rect = [
                    (j * cell_size + cell_border,
                     i * cell_size + cell_border),
                    ((j + 1) * cell_size - cell_border,
                     (i + 1) * cell_size - cell_border)
                ]
                if self.crossword.structure[i][j]:
                    draw.rectangle(rect, fill="white")
                    if letters[i][j]:
                        w, h = draw.textsize(letters[i][j], font=font)
                        draw.text(
                            (rect[0][0] + ((interior_size - w) / 2),
                             rect[0][1] + ((interior_size - h) / 2) - 10),
                            letters[i][j], fill="black", font=font
                        )

        img.save(filename)

    def solve(self):
        """
        Enforce node and arc consistency, and then solve the CSP.
        """
        self.enforce_node_consistency()
        self.ac3()
        return self.backtrack(dict())

    def enforce_node_consistency(self):
        """
        Update `self.domains` such that each variable is node-consistent.
        (Remove any values that are inconsistent with a variable's unary
         constraints; in this case, the length of the word.)
        """
        for var in self.domains.keys():
            words = self.domains[var].copy()
            for word in words:
                if len(word) != var.length:
                    self.domains[var].remove(word)

    def revise(self, x, y):
        """
        Make variable `x` arc consistent with variable `y`.
        To do so, remove values from `self.domains[x]` for which there is no
        possible corresponding value for `y` in `self.domains[y]`.

        Return True if a revision was made to the domain of `x`; return
        False if no revision was made.
        """
        revised = False
        overlap = self.crossword.overlaps[x, y]
        # if x and y are overlapping
        if overlap != None:
            # make copy of x's domain just to iterate through
            # then remove x value in its actual domain
            x_domain_copy = self.domains[x].copy()
            for x_val in x_domain_copy:
                y_domain_copy = self.domains[y].copy()
                for y_val in self.domains[y]:
                    if y_val == x_val:
                        y_domain_copy.remove(y_val)
                    else:
                        i, j = overlap
                        if x_val[i] != y_val[j]:   
                            y_domain_copy.remove(y_val)
                # if there's no value left in y's domain after enforcing arc consistent --> remove this x value
                if len(y_domain_copy) == 0:
                    self.domains[x].remove(x_val)
                    revised = True 
        return revised                            

    def ac3(self, arcs=None):
        """
        Update `self.domains` such that each variable is arc consistent.
        If `arcs` is None, begin with initial list of all arcs in the problem.
        Otherwise, use `arcs` as the initial list of arcs to make consistent.

        Return True if arc consistency is enforced and no domains are empty;
        return False if one or more domains end up empty.
        """
        # if arcs == None: queue = all arcs in problem
        if arcs == None:
            # iterate through self.crossword.overlaps.items() to get all pairs of variables in the problem with overlap value != None
            queue = [k for (k, v) in self.crossword.overlaps.items() if v != None]

        # if arcs != None: queue = arcs
        else:
            queue = arcs

        # queue - first in first out
        while len(queue) != 0:
            # dequeue queue - removing from the right
            arc = queue.pop()
            x, y = arc
            if self.revise(x, y):
                if self.domains[x] == 0:
                    return False
                x_other_neighbors = self.crossword.neighbors(x) - {y}
                for z in x_other_neighbors:
                    # enqueue - adding to the right
                    queue.append((z, x))
        return True                  

    def assignment_complete(self, assignment):
        """
        Return True if `assignment` is complete (i.e., assigns a value to each
        crossword variable); return False otherwise.
        """
        # check to see if assignment keys are the same as crossword variables
        crossword_variables = set(self.crossword.variables)
        assignment_keys = set(assignment.keys())
        if assignment_keys == crossword_variables:
            # check if all values in assignment are not empty string
            if all(value != "" for value in assignment.values()):
                return True
        return False

    def consistent(self, assignment):
        """
        Return True if `assignment` is consistent (i.e., words fit in crossword
        puzzle without conflicting characters); return False otherwise.
        """
        assignment_values = list(assignment.values())
        # if there are no duplicated values in assignment values
        if len(set(assignment_values)) == len(assignment_values):
            # if all values have correct var lengths:
            if all(len(val) == var.length for (var, val) in assignment.items()):
                # iterate through overlaps
                for key, val in self.crossword.overlaps.items():
                    v1, v2 = key
                    # if value for that pair != None and both v1 and v2 are in assignments
                    if (val != None) and (v1 in assignment.keys()) and (v2 in assignment.keys()):
                        # check if v1's index i character is the same as v2's index j character
                        i, j = val
                        if assignment[v1][i] != assignment[v2][j]:
                            return False
                return True
        return False        

    def order_domain_values(self, var, assignment):
        """
        Return a list of values in the domain of `var`, in order by
        the number of values they rule out for neighboring variables.
        The first value in the list, for example, should be the one
        that rules out the fewest values among the neighbors of `var`.
        """
        # get all of variable var's unassigned neighbors
        unassigned_variables = set(self.crossword.variables) - set(assignment.keys())
        unassigned_neighbors = [x for x in self.crossword.neighbors(var) if x in unassigned_variables]
        # create dictionary val_n with
        # key = values in var's domain
        # value = possible choices for var's unassigned neighbors (n)
        val_n = dict()
        for val in self.domains[var]:
            # initialize all n values to 0
            val_n[val] = 0
            # iterate through all unassigned neighbors
            for neighbor in unassigned_neighbors:
                # clean up neighbor's domain by removing val (if it's in there)
                neighbor_domain = set(self.domains[neighbor]) - {val}
                # make copy of neighbor's domain to loop through
                neighbor_domain_copy = neighbor_domain.copy()
                # get overlap indexes
                i, j = self.crossword.overlaps[var, neighbor]
                # remove more values from neighbor's domain 
                # if neighbor's index j character != var's index i character
                for neighbor_val in neighbor_domain_copy:
                    if neighbor_val[j] != val[i]:
                        neighbor_domain.remove(neighbor_val)
                # count the number of values left in neighbor's domain
                val_n[val] += len(neighbor_domain)
        
        # sort val_n dictionary by value in ascending order
        sorted_val_n = dict(sorted(val_n.items(), key=operator.itemgetter(1)))
        return(list(sorted_val_n.keys()))        

    def select_unassigned_variable(self, assignment):
        """
        Return an unassigned variable not already part of `assignment`.
        Choose the variable with the minimum number of remaining values
        in its domain. If there is a tie, choose the variable with the highest
        degree. If there is a tie, any of the tied variables are acceptable
        return values.
        """
        # set of unassigned variables
        unassigned_variables = set(self.crossword.variables) - set(assignment.keys())
        # create a dictionary with
        # key = unassigned variable var
        # value = [m, n] in which m = remaining values for var, and n = degrees
        var_analysis = dict()
        for var in unassigned_variables:
            var_analysis[var] = [0, 0]
            var_analysis[var][0] += len(self.domains[var])
            var_analysis[var][1] += len(self.crossword.neighbors(var))
        
        # sort var_analysis by remaining values - ascending order
        sorted_var_analysis = dict(sorted(var_analysis.items(), key=operator.itemgetter(1)))
        # get smallest remaining values in sorted val_n
        min_remaining_values = list(sorted_var_analysis.values())[0][0]
        # create sub dictionary of variables with min remaining values in domains
        min_domain = {k: v for (k, v) in sorted_var_analysis.items() if v[0] == min_remaining_values}
        # sort min_domain by degree values - ascending order
        sorted_min_domain = dict(sorted(min_domain.items(), key=operator.itemgetter(1)))
        # pick variable with highest degrees - last key
        output = list(sorted_min_domain.keys())[-1]
        
        return output

    def backtrack(self, assignment):
        """
        Using Backtracking Search, take as input a partial assignment for the
        crossword and return a complete assignment if possible to do so.

        `assignment` is a mapping from variables (keys) to words (values).

        If no assignment is possible, return None.
        """
        if self.assignment_complete(assignment):
            return assignment
        else:
            unassigned_variables = set(self.crossword.variables) - set(assignment.keys())
            var = self.select_unassigned_variable(assignment)
            queue = [(neighbor, var) for neighbor in self.crossword.neighbors(var) if not neighbor in assignment.keys()]
            # get ordered domain for this var
            ordered_domain = self.order_domain_values(var, assignment)
            # iterate through ordered domain, start with value with least-constraining value
            for value in ordered_domain:
                new_assignment = assignment.copy()
                # add {var: value} pair to assignment
                new_assignment[var] = value
                # if new assignment is consistent
                if self.consistent(new_assignment):
                    # if ac3 returns True --> arc consistency is enforced and no domains are empty
                    if self.ac3(arcs=queue):
                        # check to see if any unassigned variable has a domain with 1 value left 
                        # add those new inferences to assignment
                        for x in unassigned_variables:
                            if len(self.domains[x]) == 1 and (not x in new_assignment.keys()):
                                new_assignment[x] = list(self.domains[x])[0]
                    # run backtrack on new assignment
                    result = self.backtrack(new_assignment)
                    if result != None:
                        return result
        return None


def main():

    # Check usage
    if len(sys.argv) not in [3, 4]:
        sys.exit("Usage: python generate.py structure words [output]")

    # Parse command-line arguments
    structure = sys.argv[1]
    words = sys.argv[2]
    output = sys.argv[3] if len(sys.argv) == 4 else None

    # Generate crossword
    crossword = Crossword(structure, words)
    creator = CrosswordCreator(crossword)
    assignment = creator.solve()

    # Print result
    if assignment is None:
        print("No solution.")
    else:
        creator.print(assignment)
        if output:
            creator.save(assignment, output)


if __name__ == "__main__":
    main()
