import java.io.IOException;
import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;

import javax.swing.JFrame;
import javax.swing.JPanel;

import java.io.BufferedWriter;
import java.io.FileWriter;

import java.awt.*;


public class Ising {
    int nt      = 4;          //  number of temperature points
    int N       = 100;          //  size of the lattice, N x N
    int eqSteps = 512;        //  number of MC sweeps for equilibration
    int mcSteps = 1024;       //  number of MC sweeps for calculation
    double T1      = 0.1;
    double T2      = 0.8;
    double J = 0.5;
    double H = 0.5;

    double[] T = {0.1, 0.2, 0.3, 0.8};
    double[] E = new double[nt];
    double[] M = new double[nt];
    double[] C = new double[nt];
    double[] X = new double[nt];
    double n1 = 1.0/(mcSteps*N*N);
    double n2 = 1.0/(mcSteps*N*N);

    int[][] grid = new int[N][N];


    public Ising() {
        mainCalc();
    }

    void init() {
        Random random = new Random();
        for (int i = 0; i < this.N; i++) {
            for (int j = 0; j < this.N; j++) {
                if (Math.sqrt((N / 2 - i)*(N / 2 - i) + (N / 2 - j)*(N / 2 - j)) <= N / 6) {
                    grid[i][j] = 0;
                    continue;
                }
                int randomNumber = random.nextBoolean() ? 1 : -1;
                grid[i][j] = randomNumber;
            }
        }
    }

    void mcmove(double beta) {
        int s;
        double nb, cost;
        
        for(int i = 0; i < N; i ++) {
            for (int j = 0; j < N; j++) {
                int a = ThreadLocalRandom.current().nextInt(0, N);
                int b = ThreadLocalRandom.current().nextInt(0, N);
                s = grid[a][b];
                nb = grid[Math.floorMod((a + 1), N)][b] + grid[a][Math.floorMod((b + 1), N)] + grid[Math.floorMod((a - 1), N)][b] + grid[a][Math.floorMod((b - 1), N)];
                cost = 2 * s * nb;

                Random rand = new Random();
                if (cost < 0) {
                    s *= -1;
                }
                else if ((double)Math.abs(rand.nextInt())/Integer.MAX_VALUE < Math.exp(-cost * beta)) {
                    s *= -1;
                }

            grid[a][b] = s;
            }
        }
    }

    double calcEnergy() {
        double energy = 0.;
        for(int i = 0; i < N; i ++) {
            for (int j = 0; j < N; j++) {
                energy += -J * grid[i][j] * grid[Math.floorMod((i + 1), N)][j] - H * grid[i][j];
                energy += -J * grid[i][j] * grid[i][Math.floorMod((j + 1), N)] - H * grid[i][j];
                energy += -J * grid[i][j] * grid[Math.floorMod((i - 1), N)][j] - H * grid[i][j];
                energy += -J * grid[i][j] * grid[i][Math.floorMod((j - 1), N)] - H * grid[i][j];
                // nb = -J * (grid[(i + 1) % N][j] + grid[i][(j + 1) % N] + grid[(i - 1) % N][j] + grid[i][(j - 1) % N]);
            }
        }
        return energy / 2;
    }

    int calcMag() {
        int mag = 0;
        for(int i = 0; i < N; i ++) {
            for (int j = 0; j < N; j++) {
                mag += grid[i][j];
            }
        }
        return mag;
    }

    void mainCalc() {
        init();         // initialise
        int[][] gridCopy = new int[N][N];
        for(int i = 0; i < N; i++) 
            for(int j = 0; j < N; j++)
                gridCopy[i][j] = grid[i][j];


        for (int t = 0; t < T.length; t++) {
        
        for(int i = 0; i < N; i++) 
            for(int j = 0; j < N; j++)
                grid[i][j] = gridCopy[i][j];

        double E1 = 0; int M1 = 0; double E2 = 0; int M2 = 0;
        double iT = 1.0 / T[t];
        double iT2 = iT*iT;
        
        for (int i = 0; i < eqSteps; i++) {   // equilibrate
            mcmove(iT);                       // Monte Carlo moves
        }       
        double Ene = 0;
        int Mag = 0;
        for (int i = 0; i < mcSteps; i++) {
            mcmove(iT);         
            Ene = calcEnergy();     // calculate the energy
            Mag = calcMag();       // calculate the magnetisation
            E1 = E1 + Ene;
            M1 = M1 + Mag;
            M2 = M2 + Mag*Mag;
            E2 = E2 + Ene*Ene;
        }


        JFrame isingJFrame;
        isingJFrame = initialJFrame(grid);
        // try {
        //     Thread.sleep(5000);
        // } catch(InterruptedException e) {
        //     System.out.println("got interrupted!");
        // }

        // divide by number of sites and iteractions to obtain intensive values    
        E[t] = n1 * E1;
        M[t] = n1 * M1;
        C[t] = (n1 * E2 - n2 * E1 * E1) * iT2;
        X[t] = (n1 * M2 - n2 * M1 * M1) * iT;
        }


        try (FileWriter fileWriter = new FileWriter("E.txt")) {
            BufferedWriter outputWriter = null;
            outputWriter = new BufferedWriter(fileWriter);
            for (int i = 0; i < E.length; i++) {
                outputWriter.write(E[i]+"\n");
            }
            outputWriter.flush();  
            outputWriter.close(); 
        } catch (IOException e) {
            System.out.println(e);
        }

        try (FileWriter fileWriter = new FileWriter("M.txt")) {
            BufferedWriter outputWriter = null;
            outputWriter = new BufferedWriter(fileWriter);
            for (int i = 0; i < M.length; i++) {
                outputWriter.write(M[i]+"\n");
            }
            outputWriter.flush();  
            outputWriter.close(); 
        } catch (IOException e) {
            System.out.println(e);
        }

        try (FileWriter fileWriter = new FileWriter("C.txt")) {
            BufferedWriter outputWriter = null;
            outputWriter = new BufferedWriter(fileWriter);
            for (int i = 0; i < C.length; i++) {
                outputWriter.write(C[i]+"\n");
            }
            outputWriter.flush();  
            outputWriter.close(); 
        } catch (IOException e) {
            System.out.println(e);
        }

        try (FileWriter fileWriter = new FileWriter("X.txt")) {
            BufferedWriter outputWriter = null;
            outputWriter = new BufferedWriter(fileWriter);
            for (int i = 0; i < X.length; i++) {
                outputWriter.write(X[i]+"\n");
            }
            outputWriter.flush();  
            outputWriter.close(); 
        } catch (IOException e) {
            System.out.println(e);
        }

    }

    public JFrame initialJFrame(int[][] currentArray) {
        JFrame frame = new JFrame();
        frame.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        frame.setSize(800, 800);
        frame.setLayout(new GridLayout(this.N, this.N));
        for (int i = 0; i < this.N; i++) {
            for (int j = 0; j < this.N; j++) {
                JPanel panel = new JPanel();
                if (currentArray[i][j] == 1) {
                    panel.setBackground((new Color(100, 100, 200))); // blue
                } 
                else if (currentArray[i][j] == 0) {
                    panel.setBackground((new Color(255, 255, 255))); // white
                }
                else {
                    panel.setBackground(new Color(240, 150, 75)); // orange
                }
                frame.add(panel);
            }
        }
        frame.setVisible(true);
        return frame;
    }


}
