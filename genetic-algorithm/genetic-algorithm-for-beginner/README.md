# Genetic Algorithm for Beginner

A genetic algorithm (GA) is great for finding solutions to complex search problems. 

## How it work

Algorithm of GA:

![Process][./process.png]

Step 1: Generate any random routes and calculate their fitness
value.

Step 2: Repeat the following procedure for a given number of
iteration times:

  1. Select any two best routes from the given random routes.
  2. Reproduce those two routes to produce new best routes.
  3. After Reproduction replace the best new routes with any two worst routes.

Step 3: Return the best route.

## Pseudo Code for a Basic Genetic Algorithm

```java
generation = 0;
population[generation] = initializePopulation(populationSize);
evaluatePopulation(population[generation]);
While isTerminationConditionMet() == false do
  parents = selectParents(population[generation]);
  population[generation+1] = crossover(parents);
  population[generation+1] = mutate(population[generation+1]);
  evaluatePopulation(population[generation]);
  generation++;
End loop;
```

* The pseudo code begins with creating the genetic algorithm’s initial population.
* This population is then evaluated to find the fitness values of its individuals.
* Next, a check is run to decide if the genetic algorithm’s termination condition has been met.
* If it hasn’t, the genetic algorithm begins looping and the population goes through its first round of crossover and mutation before finally being reevaluated.
* From here, crossover and mutation are continuously applied until the termination condition is met, and the genetic algorithm terminates.


## Implement a Genetic Algorithm in Java

We will have, at minimum, four classes:

* **GeneticAlgorithm** class, which abstracts the genetic algorithm itself and provides problem-specific implementations of interface methods, such as crossover, mutation, fitness evaluation, and termination condition checking.
* **Individual** class, which represents a single candidate solution and its chromosome.
* **Population** class, which represents a population or a generation of Individuals, and applies group-level operations to them.
* **AllOnesGA** class, contains the “main” method, some bootstrap code, the concrete version of the pseudocode above, and any supporting work that a specific problem may need. These classes will be named according to the problem it solves, e.g. “AllOnesGA”, “RobotController”, etc.


###### GeneticAlgorithm.java

```java
/**
* Lots of comments in the source that are omitted here!
*/
public class GeneticAlgorithm {
  private int populationSize;
  private double mutationRate;
  private double crossoverRate;
  private int elitismCount;
  public GeneticAlgorithm(int populationSize, double mutationRate, double
    crossoverRate, int elitismCount) {
    this.populationSize = populationSize;
    this.mutationRate = mutationRate;
    this.crossoverRate = crossoverRate;
    this.elitismCount = elitismCount;
  }

  public Population initPopulation(int chromosomeLength) {
    Population population = new Population(this.populationSize, chromosomeLength);
    return population;
  }

  public double calcFitness(Individual individual) {
    // Track number of correct genes
    int correctGenes = 0;
    // Loop over individual's genes
    for (int geneIndex = 0; geneIndex < individual.getChromosomeLength(); geneIndex++) {
      // Add one fitness point for each "1" found
      if (individual.getGene(geneIndex) == 1) {
        correctGenes += 1;
      }
    }
    // Calculate fitness
    double fitness = (double) correctGenes / individual.
    getChromosomeLength();
    // Store fitness
    individual.setFitness(fitness);
    return fitness;
  }

  public void evalPopulation(Population population) {
    double populationFitness = 0;
    for (Individual individual : population.getIndividuals()) {
      populationFitness += calcFitness(individual);
    }
    population.setPopulationFitness(populationFitness);
  }

  public boolean isTerminationConditionMet(Population population) {
    for (Individual individual : population.getIndividuals()) {
      if (individual.getFitness() == 1) {
        return true;
      }
    }
    return false;
  }

  public Individual selectParent(Population population) {
    // Get individuals
    Individual individuals[] = population.getIndividuals();
    // Spin roulette wheel
    double populationFitness = population.getPopulationFitness();
    double rouletteWheelPosition = Math.random() * populationFitness;
    // Find parent
    double spinWheel = 0;
    for (Individual individual : individuals) {
      spinWheel += individual.getFitness();
      if (spinWheel >= rouletteWheelPosition) {
        return individual;
      }
    }
    return individuals[population.size() - 1];
  }

  public Population crossoverPopulation(Population population) {
    // Create new population
    Population newPopulation = new Population(population.size());
    // Loop over current population by fitness
    for (int populationIndex = 0; populationIndex < population.size(); populationIndex++) {
      Individual parent1 = population.getFittest(populationIndex);
      // Apply crossover to this individual?
      if (this.crossoverRate > Math.random() && populationIndex > this.elitismCount) {
        // Initialize offspring
        Individual offspring = new Individual(parent1.getChromosomeLength());
        // Find second parent
        Individual parent2 = selectParent(population);
        // Loop over genome
        for (int geneIndex = 0; geneIndex < parent1.getChromosomeLength(); geneIndex++) {
          // Use half of parent1's genes and half of parent2's genes
          if (0.5 > Math.random()) {
            offspring.setGene(geneIndex,
            parent1.getGene(geneIndex));
          } else {
            offspring.setGene(geneIndex,
            parent2.getGene(geneIndex));
          }
        }
        // Add offspring to new population
        newPopulation.setIndividual(populationIndex, offspring);
      } else {
        // Add individual to new population without applying crossover
        newPopulation.setIndividual(populationIndex, parent1);
      }
    }
    return newPopulation;
  }

  public Population mutatePopulation(Population population) {
    // Initialize new population
    Population newPopulation = new Population(this.populationSize);
    // Loop over current population by fitness
    for (int populationIndex = 0; populationIndex < population.size();
      populationIndex++) {
      Individual individual = population.
      getFittest(populationIndex);
      // Loop over individual's genes
      for (int geneIndex = 0; geneIndex < individual.
        getChromosomeLength(); geneIndex++) {
        // Skip mutation if this is an elite individual
        if (populationIndex >= this.elitismCount) {
          // Does this gene need mutation?
          if (this.mutationRate > Math.random()) {
            // Get new gene
            int newGene = 1;
            if (individual.getGene(geneIndex) == 1) {
              newGene = 0;
            }
          // Mutate gene
          individual.setGene(geneIndex, newGene);
        }
      }
    }
    // Add individual to population
    newPopulation.setIndividual(populationIndex, individual);
    // Return mutated population
    return newPopulation;
  }
  /**
  * Many more methods implemented later...
  */
}
```

###### Individual.java

```java
public class Individual {
  private int[] chromosome;
  private double fitness = -1;
  public Individual(int[] chromosome) {
    // Create individual chromosome
    this.chromosome = chromosome;
  }
  public Individual(int chromosomeLength) {
    this.chromosome = new int[chromosomeLength];
    for (int gene = 0; gene < chromosomeLength; gene++) {
      if (0.5 < Math.random()) {
        this.setGene(gene, 1);
      } else {
        this.setGene(gene, 0);
      }
    }
  }

  public int[] getChromosome() {
    return this.chromosome;
  }

  public int getChromosomeLength() {
    return this.chromosome.length;
  }

  public void setGene(int offset, int gene) {
    this.chromosome[offset] = gene;
  }

  public int getGene(int offset) {
    return this.chromosome[offset];
  }

  public void setFitness(double fitness) {
    this.fitness = fitness;
  }

  public double getFitness() {
    return this.fitness;
  }

  public String toString() {
    String output = "";
    for (int gene = 0; gene < this.chromosome.length; gene++) {
      output += this.chromosome[gene];
    }
    return output;
  }
}
```

###### Population.java

```java
public class Population {
  private Individual population[];
  private double populationFitness = -1;

  public Population(int populationSize) {
    this.population = new Individual[populationSize];
  }

  public Population(int populationSize, int chromosomeLength) {
    this.population = new Individual[populationSize];
    for (int individualCount = 0; individualCount <
      populationSize; individualCount++) {
      Individual individual = new
      Individual(chromosomeLength);
      this.population[individualCount] = individual;
    }
  }

  public Individual[] getIndividuals() {
    return this.population;
  }

  public Individual getFittest(int offset) {
    Arrays.sort(this.population, new Comparator<Individual>() {
      @Override
      public int compare(Individual o1, Individual o2) {
        if (o1.getFitness() > o2.getFitness()) {
          return -1;
        } else if (o1.getFitness() < o2.getFitness()) {
          return 1;
        }
        return 0;
      }
    });
    return this.population[offset];
  }

  public void setPopulationFitness(double fitness) {
    this.populationFitness = fitness;
  }

  public double getPopulationFitness() {
    return this.populationFitness;
  }

  public int size() {
    return this.population.length;
  }

  public Individual setIndividual(int offset, Individual individual) {
    return population[offset] = individual;
  }

  public Individual getIndividual(int offset) {
    return population[offset];
  }

  public void shuffle() {
    Random rnd = new Random();
    for (int i = population.length - 1; i > 0; i--) {
      int index = rnd.nextInt(i + 1);
      Individual a = population[index];
      population[index] = population[i];
      population[i] = a;
    }
  }
}
```

###### AllOnesGA.java

```java
public class AllOnesGA {
  public static void main(String[] args) {
    // Create GA object
    GeneticAlgorithm ga = new GeneticAlgorithm(100, 0.001, 0.95, 0);
    // Initialize population
    Population population = ga.initPopulation(50);
    // Evaluate population
    ga.evalPopulation(population);
    // Keep track of current generation
    int generation = 1;
    while (ga.isTerminationConditionMet(population) == false) {
      // Print fittest individual from population
      System.out.println("Best solution: " + population.
      getFittest(0).toString());
      // Apply crossover
      population = ga.crossoverPopulation(population);
      // Apply mutation
      population = ga.mutatePopulation(population);
      // Evaluate population
      ga.evalPopulation(population);
      // Increment the current generation
      generation++;
    }
    System.out.println("Found solution in " + generation + "generations");
    System.out.println("Best solution: " + population.getFittest(0).toString());
  }
}
```

Output may look like:

```
Best solution: 11001110100110111111010111001001100111110011111111
Best solution: 11001110100110111111010111001001100111110011111111
Best solution: 11001110100110111111010111001001100111110011111111
[ ... Lots of lines omitted here ... ]
Best solution: 11111111111111111111111111111011111111111111111111
Best solution: 11111111111111111111111111111011111111111111111111
Found solution in 113 generations
Best solution: 11111111111111111111111111111111111111111111111111
```
