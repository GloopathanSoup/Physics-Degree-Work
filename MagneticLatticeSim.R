#MAGNETIC LATTICE SIM

#This was a part of one of my second year projects, creating a simplified simulation of a 2D magnetic lattice
#Officially it received a perfect score, but in retrospect I was lucky that they didn't penalise inefficiency and excessive commenting
#The implementation is simpler than it may sound to a non-physicist
#It starts off as a random 2D array of 1s and -1s
#It progresses according to simple physical rules relating to entropy, with spin flips (*-1) occuring randomly, but being more likely when the flip would result in a lower energy (more aligned) state for the whole lattice
#From this simple model, interesting phenomena emerge, mirroring real physics
#At low input temperatures (temper < 2.5), it tends toward larger aligned areas (eventually the whole aray aligning to 1 or -1), ordered like the structure of a solid
#At high input temperatures (temper > 2.5), it remains in total random chaos without structure, like a gas


N.steps <- function(N, temper = 1, n = 50, boltz = 1){
  
  S.dim <- c(n, n) #defining the dimensions of the lattice for use in the array() function
  
  S <- array(1, dim = S.dim) #defines the function as a 2 dimensional array of dimensions defined above with sides n (where 1 is a placeholder value)
  
  pos.spins <- c(-1, 1) #defines possible spin values
  
  initial.spins <- sample(pos.spins, size=n^2, replace = TRUE) #creates an array defining the initial spin values as random replaced samples of the possible spins (-1 and 1) of a size that of the size of the lattice (n^2)
  
  S[1:n, 1:n] <- initial.spins #applies the initial spins array onto the n*n lattice array
  
  #defining i and j as the range of non-edge elements of the lattice
  i <- (2:(n-1))
  j <- (2:(n-1))
  
  energy <- function(S, J = 1) { #creates a function outputting the total energy over the lattice S for a given minimal unit of energy J with a default value of 1
    
    #defining the energy of a single non-edge element
    NE.element.energy <- J*((-1)*S[i, j]*S[(i-1), j] + (-1)*S[i, j]*S[(i+1), j] + (-1)*S[i, j]*S[i, (j-1)] + (-1)*S[i, j]*S[i, (j+1)]) #energy of s[i, j] computed as the sum of energy as a result of all adjacent elements
    
    sum.NE.edge.energy <- sum(NE.element.energy) #sums the element energy over every non-edge element
    
    #defining the energy over the -i edge
    ni.element.energy <- J*((-1)*S[1, j]*S[2, j] + (-1)*S[1, j]*S[n, j] #the parts of energy change defined by the -i edge's adjacency to the i=2 row and the opposing i=n row 
                            + (-1)*S[1, j]*S[1, (j+1)] + (-1)*S[1, j]*S[1, (j-1)]) #the parts of energy change defined by elements in the -i edge and true adjacent elements in the j and -j directions
    
    #defining the energy over the +i edge
    pi.element.energy <- J*((-1)*S[n, j]*S[(n-1), j] + (-1)*S[n, j]*S[1, j] #the parts of energy change defined by the +i edge's adjacency to the i=(n-1) row and the opposing i=1 row 
                            + (-1)*S[n, j]*S[n, (j+1)] + (-1)*S[n, j]*S[n, (j-1)]) #the parts of energy change defined by elements in the i edge and true adjacent elements in the j and -j directions
    
    #defining the energy over the -j edge
    nj.element.energy <- J*((-1)*S[i, 1]*S[i, 2] + (-1)*S[i, 1]*S[i, n] #the parts of energy change defined by the -j edge's adjacency to the j=2 row and the opposing j=n row 
                            + (-1)*S[i, 1]*S[(i+1), 1] + (-1)*S[i, 1]*S[(i-1), 1]) #the parts of energy change defined by elements in the -j edge and true adjacent elements in the i and -i directions
    
    #defining the energy over the +j edge
    pj.element.energy <- J*((-1)*S[i, n]*S[i, (n-1)] + (-1)*S[i, n]*S[i, 1] #the parts of energy change defined by the +j edge's adjacency to the i=(n-1) row and the opposing i=1 row 
                            + (-1)*S[i, n]*S[(i+1), n] + (-1)*S[i, n]*S[(i-1), n]) #the parts of energy change defined by elements in the j edge and true adjacent elements in the i and -i directions
    
    #defining the energy of the [1, 1] corner element
    oo.element.energy <- J*((-1)*S[1, 1]*S[1, n] + (-1)*S[1, 1]*S[n, 1] #the parts of energy change defined by S[1, 1] and the opposite corners S[1, n] and S[n, 1]
                            + (-1)*S[1, 1]*S[1, 2] + (-1)*S[1, 1]*S[2, 1]) #the parts of energy change defined by S[1, 1] and its true adjacent elements
    
    #defining the energy of the [n, 1] corner element
    no.element.energy <- J*((-1)*S[n, 1]*S[n, n] + (-1)*S[n, 1]*S[1, 1] #the parts of energy change defined by S[n, 1] and the opposite corners S[n, n] and S[1, 1]
                            + (-1)*S[n, 1]*S[n, 2] + (-1)*S[n, 1]*S[n, (n-1)]) #the parts of energy change defined by S[n, 1] and its true adjacent elements
    
    #defining the energy of the [1, n] corner element
    on.element.energy <- J*((-1)*S[1, n]*S[1, 1] + (-1)*S[1, n]*S[n, n] #the parts of energy change defined by S[1, n] and the opposite corners S[1, 1] and S[n, n]
                            + (-1)*S[1, n]*S[1, (n-1)] + (-1)*S[1, n]*S[(n-1), 1]) #the parts of energy change defined by S[1, n] and its true adjacent elements
    
    #defining the energy of the [n, n] corner element
    nn.element.energy <- J*((-1)*S[n, n]*S[1, n] + (-1)*S[n, n]*S[n, 1] #the parts of energy change defined by S[n, n] and the opposite corners S[1, n] and S[n, 1]
                            + (-1)*S[n, n]*S[n, (n-1)] + (-1)*S[n, n]*S[(n-1), n]) #the parts of energy change defined by S[n, n] and its true adjacent elements
    
    E <- (sum(NE.element.energy) + sum(ni.element.energy) + sum(pi.element.energy) + sum(nj.element.energy) + sum(pj.element.energy) + oo.element.energy + on.element.energy + no.element.energy + nn.element.energy)/2 #sums energy of every element into one object and divides by 2 to account for each pair being computed twice
    
    return(E)
    
  }
  
  
  #a function to compute the energy change to the system given the spin flip of one element
  dif.E <- function(k, l, J = 1) { #sets the function as depending on the position of the flipped element [k, l] and on the energy constant, set to a default of 1 as before
    
    if(1 < k && k < n && 1 < l && l < n){ #sets the following code to be computed only if a non-edge element is inputted (if k and l are between 1 and n, exclusive)
      dE <- J*(S[k, l]*S[(k-1), l] + S[k, l]*S[(k+1), l] + S[k, l]*S[k, (l-1)] + S[k, l]*S[k, (l+1)]) #energy of S[k, l] after a spin flip computed as the sum of energy as a result of all adjacent elements
      - J*((-1)*S[k, l]*S[(k-1), l] + (-1)*S[k, l]*S[(k+1), l] + (-1)*S[k, l]*S[k, (l-1)] + (-1)*S[k, l]*S[k, (l+1)])} #energy of s[k, l] before the flip, subtracted from the energy after to give the change in energy
    
    if(k == 1 && l != 1 && l != n){ #sets the following code to be computed only if a -k edge element is inputted (if k = 1 and isn't a corner element of that edge)
      dE <- J*(S[1, l]*S[2, l] + S[1, l]*S[n, l] + S[1, l]*S[1, (l+1)] + S[1, l]*S[1, (l-1)]) #energy after the flip
      - J*((-1)*S[1, l]*S[2, l] + (-1)*S[1, l]*S[n, l] + (-1)*S[1, l]*S[1, (l+1)] + (-1)*S[1, l]*S[1, (l-1)])} #energy before the flip, subtracted from the energy after to give the change in energy
    
    if(k == n && l != 1 && l != n){ #sets the following code to be computed only if a +k edge element is inputted
      dE <- J*(S[n, l]*S[(n-1), l] + S[n, l]*S[1, l] + S[n, l]*S[n, (l+1)] + S[n, l]*S[n, (l-1)]) #energy after the flip
      - J*((-1)*S[n, l]*S[(n-1), l] + (-1)*S[n, l]*S[1, l] + (-1)*S[n, l]*S[n, (l+1)] + (-1)*S[n, l]*S[n, (l-1)])} #energy before the flip, subtracted from the energy after to give the change in energy
    
    if(l == 1 && k != 1 && k != n){ #sets the following code to be computed only if a -l edge element is inputted
      dE <- J*(S[k, 1]*S[k, 2] + S[k, 2]*S[k, n] + S[k, 1]*S[(k+1), 1] + S[k, 1]*S[(k-1), 1]) #energy after the flip
      - J*((-1)*S[k, 1]*S[k, 2] + (-1)*S[k, 1]*S[k, n] + (-1)*S[k, 1]*S[(k+1), 1] + (-1)*S[k, 1]*S[(k-1), 1])} #energy before the flip, subtracted from the energy after to give the change in energy
    
    if(l == n && k != 1 && k != n){ #sets the following code to be computed only if a +l edge element is inputted
      dE <- J*(S[k, n]*S[k, (n-1)] + S[k, n]*S[k, 1] + S[k, n]*S[(k+1), n] + S[k, n]*S[(k-1), n]) #energy after the flip
      - J*((-1)*S[k, n]*S[k, (n-1)] + (-1)*S[k, n]*S[k, 1] + (-1)*S[k, n]*S[(k+1), n] + (-1)*S[k, n]*S[(k-1), n])} #energy before the flip, subtracted from the energy after to give the change in energy
    
    if(k == 1 && l == 1){ #sets the following code to be computed only if the [1, 1] element is inputted
      dE <- J*(S[n, 1]*S[n, n] + S[n, 1]*S[1, 1] + S[n, 1]*S[n, 2] + S[n, 1]*S[n, (n-1)]) #energy after the flip
      - J*((-1)*S[n, 1]*S[n, n] + (-1)*S[n, 1]*S[1, 1] + (-1)*S[n, 1]*S[n, 2] + (-1)*S[n, 1]*S[n, (n-1)])}#energy before the flip, subtracted from the energy after to give the change in energy
    
    if(k == n && l == 1){ #sets the following code to be computed only if the [n, 1] element is inputted
      dE <- J*(S[n, 1]*S[n, n] + S[n, 1]*S[1, 1] + S[n, 1]*S[n, 2] + S[n, 1]*S[n, (n-1)]) #energy after the flip
      - J*((-1)*S[n, 1]*S[n, n] + (-1)*S[n, 1]*S[1, 1] + (-1)*S[n, 1]*S[n, 2] + (-1)*S[n, 1]*S[n, (n-1)])}#energy before the flip, subtracted from the energy after to give the change in energy
    
    if(k == 1 && l == n){ #sets the following code to be computed only if the [1, n] element is inputted
      dE <- J*(S[1, n]*S[1, 1] + S[1, n]*S[n, n] + S[1, n]*S[1, (n-1)] + S[1, n]*S[(n-1), 1]) #energy after the flip
      - J*((-1)*S[1, n]*S[1, 1] + (-1)*S[1, n]*S[n, n] + (-1)*S[1, n]*S[1, (n-1)] + (-1)*S[1, n]*S[(n-1), 1])}#energy before the flip, subtracted from the energy after to give the change in energy
    
    if(k == n && l == n){ #sets the following code to be computed only if the [n, n] element is inputted
      dE <- J*(S[n, n]*S[1, n] + S[n, n]*S[n, 1] + S[n, n]*S[n, (n-1)] + S[n, n]*S[(n-1), n]) #energy after the flip
      - J*((-1)*S[n, n]*S[1, n] + (-1)*S[n, n]*S[n, 1] + (-1)*S[n, n]*S[n, (n-1)] + (-1)*S[n, n]*S[(n-1), n])}#energy before the flip, subtracted from the energy after to give the change in energy
    
    return(2*dE) #returns the difference in energy with an arbitrary 2* value which is necessary to reach expected values for unknown reasons
    
  }
  
  beta <- 1/(boltz*temper) #defines beta
  
  #placeholder arrays to store magnetisation and average energy at every step of the function
  m.list <- (1:(N + 1))
  e.list <- (1:(N + 1))
  
  t <- (1:(N + 1)) #creates values for time at each step so that energy and magnetisation can be plotted against time
  
  #records the magnetisation and average energy of the initial condition of the lattice in the first position of the list
  m.list[1] <- mean(S)
  e.list[1] <- energy(S)/(n^2)
  
  for(m in 1:N){ #loops the code of the function N times, choosing N random elements and flipping or not flipping them if appropriate
    
    i.r <- sample.int(n, size = 1, replace = TRUE) #randomises the i position of the chosen S element
    j.r <- sample.int(n, size = 1, replace = TRUE) #randomises the j position of the chosen S element
    dE.step <- dif.E(i.r, j.r)
    
    accept <- FALSE #sets the flip to be rejected if conditions aren't met
    change.energy.step <- 0 #defines the change of energy as a result of no flip occurring
    
    if(dE.step < 0){ #sets the following code to run if the flip of the element would decrease the energy of the lattice
      accept <- TRUE #sets the flip to be accepted if it would decrease the energy of the system
      S[i.r, j.r] <- (-1)*S[i.r, j.r] #flips the random element by multiplying it by -1
    } else { #sets the following code to run if the flip doesn't decrease the energy of the lattice
      u <- runif(1) #generates a random number between 0 and 1 for use in a probability distribution
      w <- exp(-(beta*dE.step)) #sets an object exponentially dependent on the energy change and beta for use in an exponential probability distribution
      if(w >= u){ #combines u and w into the previously described probability distribution
        accept <- TRUE #sets the flip to be randomly accepted based off of the previously described distribution
        S[i.r, j.r] <- (-1)*S[i.r, j.r] #flips the random element by multiplying it by -1
      }
      
      
    }
    #inserts values into the m and e lists on each instance of the loop into the correct position (m + 1 to account for the first state)
    m.list[m + 1] <- mean(S)
    if(accept == TRUE){ e.list[m + 1] <- dif.E(i.r, j.r) } else { e.list[m + 1] <- 0 }
  }
  
  #graphical outputs
  plot(t, m.list, xlab = "Time step", ylab = "Lattice magnetisation", type = "l")
  plot(t, cumsum(e.list), xlab = "Time step", ylab = "Lattice energy", type = "l")
  image(S)
  
  if(N < 50){
    #text output
    return(cat("Final state of the lattice \n", S, "\n \n Average energy at each successive state \n", cumsum(e.list),"\n \n Lattice magnetisation at each successive state \n", m.list))
  }
}
