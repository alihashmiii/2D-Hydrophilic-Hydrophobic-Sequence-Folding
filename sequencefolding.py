import itertools, random, pylab, matplotlib.pyplot as plt, math, os     # importing the essential libraries

print "the code can be used to find the optimal (thermodynamically favourable) structures of a Hydrophobic-Polar amino-acid sequence,\
confined to a 2-D lattice using a brute force approach. Please change the energy threshold on line 66 as desired", '\n'

output_dir = 'HP_sequence_to_struct'                            # saves all the images in the directory: 'HP_sequence_to_struct'
if not os.path.exists(output_dir): os.makedirs(output_dir)

base_num = raw_input('enter the number of bases: ')             # prompts the user for the number of bases forming the protein HP sequence
base_num = int(base_num)

sequence = raw_input('enter the sequence of %d HP: ' %base_num) # prompts the user for the sequence of choice as H and P's
sequence = list(sequence)                                       # breaks HPPPHPPP into ['H','P','P','P','H','P','P','P'] (string to characters)

color_seq = ['r' if x == 'P' else 'b' for x in sequence]        # assigns color to the result of the previous line (blue for H and red for P)

# making the grid (later used for plotting) and a record keeper for the places taken in the grid

grid_siz_x, grid_siz_y = base_num*2 + 2, base_num*2 + 2         # determines the size of the grid based on the number of bases entered earlier
taken = {x:0 for x in range(1,grid_siz_x*grid_siz_y + 1)}                                               # creates a #map with keys equal to the size of the grid
coordinate = [(x,y) for y in [i for i in xrange(grid_siz_x)] for x in [i for i in xrange(grid_siz_y)]]  # creates coordinates that can be assigned to available positions
keys = range(1,grid_siz_x*grid_siz_y + 1)                                                               # generates a set of all keys
grid = {}                                       # creates an empty grid #map

a = itertools.starmap(lambda x,(w,y): grid.update({x:(w,y)}), zip(keys,coordinate))         # updates the #map with keys assigned to coordinates i.e. key:(x,y)
a = list(a)

# setting the moves for moving from one grid point to the adjacent neighbours

up, down, right, left = grid_siz_x,-grid_siz_x,1,-1
move = [up,down,right,left]

# fixing the 1st and the 2nd base and also ensuring that the #map 'taken' is set to 1 for the indices of the two bases

start_pt_ind = (grid_siz_x*grid_siz_y)/2 - (grid_siz_x/2)   # setting the position of the first base
start_point_coord = grid[start_pt_ind]                      # storing coordinates of the 1st base

taken[start_pt_ind] = 1                             # assigning value of key as 1 in 'taken' after setting the first base

second_pt_ind = start_pt_ind + random.choice(move)  # randomly setting the position of the second base in the grid
second_pt_coord = grid[second_pt_ind ]              # recording the coordinate of the second base

taken[second_pt_ind] = 1                            # assigning value of key as 1 in 'taken' after setting the second base

pos = [grid[start_pt_ind],grid[second_pt_ind]]      # stores the coordinates of the first and second base

taken_init = taken.copy()                           # creates a Copy of the #map 'taken' after the indices (locations in the grid) for the first and second base are set to 1


# turn the remaining string of (base_num - 3) HP's to get all possible permutations of HPs

combinations = 4**base_num      # all possible permutations of the sequence (each H or P in the sequence has the tendency to take one of 4 positions: up, down, right, left)
store_quart_list = []           # stores the positions of all the bases in the empty list such as [0,2,3,2,1, ... 0,2,2,3]

base = base_num - 3             # number of bases that remain after the first three bases are fixed

def ConvertIntToQuart(j,base):  # the function computes all possible moves that the remaining bases can make in terms of 0, 1, 2, 3
    L = [0]*base                # creates a list with a size equal to remaining bases and values set to 0
    remainder = j

    for n in range(base-1,0,-1):                        # computing 0,1,2,3 for each base in the sequence. note: (total bases - 3)
        L[n-1] = math.floor(remainder/(4**(n-1)))
        remainder = remainder - L[n-1]*(4**(n-1))
    return L


for k in xrange(combinations):                              # calls the function ConvertIntToQuart and the result is appended to a list, store_quart_list
    store_quart_list.append(ConvertIntToQuart(k,base))

# <<<<<<<< THRESHOLD >>>>>>>>>>   set as desired

threshold = 4

# fixing the coordinates of the third base

for i in range(1,3):                    # the third base has two alternative routes to take to avoid a cross-over so we need to loop the entire structure twice

    save_val = []
    third_pt_coord = []                 # stores the coordinates of the 3rd base

    if i == 1:  # first structural iteration executed on the first run
        if (pos[0][0] == pos[1][0] and pos[1][1]>pos[0][1]) or (pos[0][1] == pos[1][1] and pos[1][0]>pos[0][0]):    # checks the position of the first two indices and assigns an index to the 3rd base
            ind = second_pt_ind + right

        elif (pos[0][0] == pos[1][0] and pos[1][1]<pos[0][1]) or (pos[0][1] == pos[1][1] and pos[1][0]<pos[0][0]):
            ind = second_pt_ind + left

    if i == 2: # second structural iteration executed on the second run
        if (pos[0][0] == pos[1][0] and pos[1][1]>pos[0][1]) or (pos[0][1] == pos[1][1] and pos[1][0]>pos[0][0]):    # checks the position of the first two indices and assigns an index to the 3rd base
            ind = second_pt_ind + up

        elif (pos[0][0] == pos[1][0] and pos[1][1]<pos[0][1]) or (pos[0][1] == pos[1][1] and pos[1][0]<pos[0][0]):
            ind = second_pt_ind + down

    third_pt_coord = grid[ind]      # coordinates of the third base are recorded
    third_pt_ind = ind

    #print i, 'taken', taken, '\n', 'taken_init', taken_init, '\n'

    # looping through the possible permutations of the arrangements that the (n-3) bases can assume
    # and passing an arrangement into the grid ("taken")

    j = 0   # initializes a counter for the length of the store_quart_list (which contains all the permutations of the positions of the bases)
    while j < len(store_quart_list):
        j+=1                            # increments the counter
        ind = third_pt_ind              # ind gets updated to the index of the third node (which can have two positions)
        taken = taken_init.copy()       # makes a copy of the grid that is used to keep track of the indices that are 'taken' or have a value 1
        taken[third_pt_ind] = 1         # updates the value of 'taken' at ind (index of 3rd node) as 1

        move_bases = store_quart_list[j-1]      # takes one of the many arrangements of HPHH... in space (a sequence of 0,1,2,3 representing the position of the node relative to its
                                                # adjacent neighbour) from store_quart_list

        old_ind = ind
        ind_save = [start_pt_ind,second_pt_ind,third_pt_ind]        # saves the first three indices into a list 'ind_save'

        ind_polarity_save = {}                              # creates an empty dictionary 'ind_polarity_save'
        pos_ind_save = {}                                   # creates an empty dictionary 'pos_ind_save'

        for k in move_bases:                # can be any element of (0,1,2,3)
            old_ind = ind                    # stores the old value of ind (checks whether ind changed or not)
            if k == 0 and taken[ind+up] == 0:   # checks if k is 0 && whether there is a node already present above the current node (false if taken at ind+up is 0)
                ind = ind + up                      # changes the value of ind to be ind + up
                taken[ind] = 1                      # changes the value of taken at ind to 1 if new node is not already taken
            elif k == 1 and taken[ind+right] == 0:      # checks if k is 1 && whether there is a node already present above the current node (false if taken at ind+right is 0)
                ind = ind + right                           # changes the value of ind to be ind + right
                taken[ind] = 1                              # changes the value of taken at ind to 1 if new node is not already taken
            elif k == 2 and taken[ind+down] == 0:   # checks if k is 2 && whether there is a node already present above the current node (false if taken at ind+down is 0)
                ind = ind + down                            # changes the value of ind to be ind + down
                taken[ind] = 1                              # changes the value of taken at ind to 1 if new node is not already taken
            elif k == 3 and taken[ind+left] == 0:   # checks if k is 3 && whether there is a node already present above the current node (false if taken at ind+left is 0)
                ind = ind + left                    # changes the value of ind to be ind + left
                taken[ind] = 1                      # changes the value of taken at ind to 1 if new node is not already taken

            ind_save.append(ind)                # saves the 4th, 5th index ... in 'ind_save'

            if ind == old_ind:                  # clears 'ind_save' if the index of the new node is already occupied and breaks the for-loop
                ind_save = []
                break
        if ind_save == []:                      # if the ind_save is empty then no need to execute the remaining portion of the while loop. instead continue to a new sequence
            continue                            # from store_quart_list

        x = [n for n in range(len(sequence))]   # gives a list of nodes [1,2, ... len(sequence)]
        a = itertools.starmap(lambda x,w,y: ind_polarity_save.update({x:(w,y)}), zip(x,ind_save,sequence))  # update ind_polarity_save
        list(a)
        b = itertools.starmap(lambda x,w,y: pos_ind_save.update({x:(w,y)}), zip(ind_save,x,sequence))       # update pos_ind_save
        list(b)

        # <<<<<<<< computing neighbour and interaction matrices for subsequent energy calculations >>>>>>>>

        neighbours = []     # stores a list of neighbour which itself is a list of four members
        for l in range(len(ind_polarity_save)):     # goes over the length of the sequence
            if ind_polarity_save[l][1] == 'H':      # checks and executes the block of code only if the node is hydrophobic
                neighbour = [0]*4                   # creates an empty list of four neighbours (since every neighbour has four adjacent neighbours)
                index = ind_polarity_save[l][0]     # index stores the index of the l (L) th node in ind_polarity_save
                for n in range(4):                  # to go over the four neighbours
                    new_index = index + move[n]     # (index + up, index + down, index + right, index + left) yields new indices
                    if taken[new_index] == 1:       # if neighbour of a particular node does exit i.e. taken[index of the neighbour] = 1
                        neighbour[n] = new_index    # replaces the member of neighbour with the index of the neighbour if neighbour exists
                    elif taken[new_index] == 0:     # if neighbour does not exit i.e. taken[index of the neighbour] = 0
                        neighbour[n] = 0            # puts 0 in the 'neighbour (list)' if the neighbour of a particular node does not exist
                neighbours.append(neighbour)        # stores the list of four neighbours of each node into a list i.e. [[neigh1,neigh2,neigh3,neigh4]_node1 , [neigh1,neigh2,neigh3,neigh4]_node2 ...]

        interaction = []                            # stores a list of base_interaction, which itself is a list of four members each delineating the interaction of the node with its neighbour
        for l in range(len(ind_polarity_save)):     # goes over the length of teh sequence
            if ind_polarity_save[l][1] == 'H':      # checks and executes the block of code only if the node is hydrophobic
                index = ind_polarity_save[l][0]     # index stores the index of the l (L) th node in ind_polarity_save
                base_interaction = [0]*4            # creates an empty list of four interactions (since every neighbour has four adjacent neighbours)
                for n in range(4):                  # to go over the four neighbours
                    neigh_index = index + move[n]   # (index + up, index + down, index + right, index + left) yields new indices
                    if taken[neigh_index] == 0:     # if neighbour of a particular node does exit i.e. taken[index of the neighbour] = 1
                        base_interaction[n] = 1     # replaces the member of base_interaction with 1 if interaction exists
                    elif pos_ind_save[neigh_index][1] == 'H' and abs(pos_ind_save[neigh_index][0] - pos_ind_save[index][0]) > 1:
                        base_interaction[n] = 1
                interaction.append(base_interaction)  # stores the list of four interactions of each node in a list i.e. [[inter1,inter2,inter3,inter4]_node1 , [inter1,inter2,inter3,inter4]_node2 ...]

        hyd_int_energy = 0                          # the initial energy of hydration (between two hydrophobic regions)
        for l in range(len(neighbours)):            # looping over the entire sequence
            hyd_count = 0                           # initial hydrophobic-hydrophobic interactions
            for n in range(4):                      # looping over the four neighbours
                if neighbours[l][n] != 0 and interaction[l][n] != 0:    # checking if the neighbours are non-zero i.e. the neighbour exists and the neighbour is a hydrophobic node
                    hyd_count+=1
            hyd_int_energy += -2*hyd_count                              # calculating the total energy due to H-H interacting pairs

        total_energy = hyd_int_energy/2         # total interaction between the hyrophobic-hydrophobic interactions

        solvent_int = [neighbours[o].count(0) for o in range(len(neighbours))]  # counts the total number of solvent interactions for each node
        num_solvent_int = sum([1 if x>1 else x for x in solvent_int])   # calculates the sum of the total interactions due to solvent for all nodes

        total_energy = total_energy + num_solvent_int       # calculates the total energy of the system

        if total_energy >= threshold:   # do not save if the total energy of the structure is greater than the threshold
            continue

        #print j,'\n','total energy', total_energy, '\n', ind_polarity_save, '\n', 'interaction', interaction, '\n', 'neighbours', neighbours,'\n', taken, '\n'

        if total_energy < threshold:            # save the structure if the total energy is less than the threshold
            save_val.append([[j,total_energy],ind_polarity_save])


    # <<< PLOTTING SEQUENCES >>>

    sigma = 0.2     # radius of the circle (represents a base)
    for k in range(len(save_val)):
        fig = plt.figure()              # creates 'k' plots for the various sequences
        pylab.axis([0, grid_siz_x, 0, grid_siz_y])      # demarcates axis for the plots

        indices = [save_val[k][1][p][0] for p in range(len(save_val[0][1]))]    # gets indices for all the nodes for a given sequence
        polarities = [save_val[k][1][p][1] for p in range(len(save_val[0][1]))] # gets polarities for all the nodes i.e. 'H' and 'P' for a particular sequence
        colour = ['r' if x == 'P' else 'b' for x in polarities]                 # assigns colours to the hydrophobic and hydrophilic regions

        pos = [grid[x] for x in indices]                                        # gets (x,y) coordinates for all the nodes using the indices

        n = 0       # local counter for while loop
        while n < base_num - 1:                 # creates the lines (edges) connecting the nodes ('H' and 'P' regions)
            if pos[n][0] == pos[n+1][0]: # draw verical lines (same x-coordinates)
                pylab.plot([pos[n][0],pos[n+1][0]],[min(pos[n][1],pos[n+1][1]),max(pos[n][1],pos[n+1][1])],'k')
            elif pos[n][1] == pos[n+1][1]: # draw horizontal lines (same y-coordinates)
                pylab.plot([min(pos[n][0],pos[n+1][0]),max(pos[n][0],pos[n+1][0])],[pos[n][1],pos[n+1][1]],'k')
            n+=1    # increments counter for the while loop

        for (x,y),c in zip(pos,colour):                         # zips the coordinates and colours of the 'H' and 'P' regions together
            circle = pylab.Circle((x, y), radius=sigma, fc=c)   # creates a circle at the given coordinates with the given radius and colour
            pylab.gca().add_patch(circle)                       # adds a circle to the same plot using gca() (getcurrentaxis)
            plt.draw()                                          # draws the circle

        if i == 1:                                              # puts the label on the graphs for the first run when the nodes are say o---o---o horizontal
            track = len(save_val)
            plt.title('total energy of structure: %d (AU)' %save_val[k][0][1])

        if i == 2:                                              # puts the label on the plots for the second run when nodes are say arranged as follows: o--o
            k += track                                                                                                                                   #  |
            plt.title('total energy of structure: %d (AU)' %save_val[k-track][0][1])                                                                     #  o

        pylab.savefig(os.path.join(output_dir, '%d.png' %k), transparent = True)        # saves the figures to a separate folder in the directory where the .py file is located
        plt.close()     # closes the plot
