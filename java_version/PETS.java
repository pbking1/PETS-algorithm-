import Jama.Matrix;
import java.lang.Math;
import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.text.*;
import java.io.BufferedReader;
import java.io.FileReader;

public class PETS{
    public static ArrayList<String> proteinName = new ArrayList<String>();
    public static ArrayList<String> drugName = new ArrayList<String>();
    public static double [][]adj;
    public static Matrix adj_M;
	public static double [][]drugVector;
	public static Matrix drugVector_M;
    public static double [][]Expression;
	public static Matrix Expression_M;
    public static double [][]Rp;
	public static Matrix Rp_M;
  	public static double boostingFactor = 1.5;
	public static double sigma = 0.8;
	public static int maxIteration = 150;
	public static double PETScoreMatrix;
	public static Matrix scoreArray_M;
	public static double [][]sMat;
	public static int numOfPro = 242;
	public static double [][]target;
	public static Matrix target_M;

	public static void InitDoubleArray(double [][]t){
		for(int i = 0; i < t.length; i++){
			for(int j = 0; j < t[i].length; j++){
				t[i][j] = 0.00;
			}
		}
	}

	public static Matrix CreateRandomMatrix(int n){
		double [][]M = new double[n][n];
		for(int i = 0; i < n; i++){
			for(int j = 0; j < n; j++){
				M[i][j] = i + j;
			}
		}
		return new Matrix(M);
	}

    public static Matrix CreateEyeMatrix(int n){
        double [][] M = new double[n][n];
        for(int i = 0; i < n; i++){
            for(int j = 0; j < n; j++){
                if(i == j)
                    M[i][j] = 1;
                else
                    M[i][j] = 0;
            }
        }
        return new Matrix(M);
    }

	public static Matrix CreateZeroMatrix(int n, int m){
		Matrix M2 = new Matrix(n, m, 0.00);
		return M2;
	}

	public static Matrix CreateSignMatrix(Matrix M3){
		double [][] M_temp = M3.getArray();
		for(int i = 0; i < M_temp.length; i++){
			for(int j = 0; j < M_temp[i].length; j++){
				if(M_temp[i][j] > 0){
					M_temp[i][j] = 1;
				}else if(M_temp[i][j] == 0){
					M_temp[i][j] = 0;
				}else if(M_temp[i][j] < 0){
					M_temp[i][j] = -1;
				}
			}
		}
		Matrix result = new Matrix(M_temp);
		return result;
	}
    
    public static Matrix AbsMatrix(Matrix M){
        double [][] M_temp = M.getArray();
        for(int i = 0; i < M_temp.length; i++){
            for(int j = 0; j < M_temp[i].length; j++){
                M_temp[i][j] = Math.abs(M_temp[i][j]);
            }
        }
        Matrix result = new Matrix(M_temp);
        return result;
    }
    
    public static Matrix DivideMatrix(Matrix M, double n){
        double [][] M_temp = M.getArray();
        for(int i = 0; i < M_temp.length; i++){
            for(int j = 0; j < M_temp[i].length; j++){
                M_temp[i][j] = M_temp[i][j] / n;
            }
        }
        Matrix result = new Matrix(M_temp);
        return result;
    }
    
    public static void PrintMatrix(Matrix M){
        double [][] M_1 = M.getArray();
        for(int i = 0; i < M_1.length; i++){
            for(int j = 0; j < M_1[i].length; j++){
                System.out.print(M_1[i][j] + " ");
            }
            System.out.print("\n");
        }
    }
    
    public static void PrintNameList(ArrayList<String> name){
        for(int i = 0; i < name.size(); i++){
            System.out.println(name.get(i));
       	}
    }
    
    public static void ReadNameFile(String filename, ArrayList<String> result) {
        try{
            BufferedReader reader = new BufferedReader(new FileReader(filename));
            String line = null;
            int j = 0;
            while((line = reader.readLine()) != null){
                //System.out.println(line);
                result.add(j, line);
                j++;
            }
        }catch(Exception e){
            e.printStackTrace();;
        }
        
    }
    
    public static double[][] ReadMatrixCsv(String name, int row, int col){
        double [][]M1 = new double[row][col];
        //System.out.println(name);
        try{
            BufferedReader reader = new BufferedReader(new FileReader(name));
            String line = null;
            int j = 0;
            while((line = reader.readLine()) != null){
                String item[] = line.split(",");
                //System.out.println(item.length);
                for(int i = 0; i < item.length; i++){
                    String temp = item[i];
                    //System.out.print(temp + " ");
                    //System.out.println(temp.length() - temp.indexOf("."));
                    if(temp.length() - temp.indexOf(".") > 2){
                        String result = temp.substring(0, temp.indexOf(".") + 3);
                        M1[j][i] = Double.parseDouble(result);
                        //System.out.print(result);
                        //System.out.print(" ");
                    }else if(temp.length() - temp.indexOf(".") < 2 && temp.indexOf("-") != 0 && temp.equals("0") == false){
                        String result = temp.substring(0, temp.indexOf(".") + 1);
                        M1[j][i] = Double.parseDouble(result);
                        //System.out.print(result);
                        //System.out.print(" ");
                    }else if(temp.equals("0") == true || temp.equals("0.00") == true || temp.indexOf("-") == 0){
                        String result = temp;
                        M1[j][i] = Double.parseDouble(result);
                        //System.out.print(result);
                        //System.out.print(" ");
                    }else{
                        String result = temp;
                        M1[j][i] = Double.parseDouble(result);
                        //System.out.print(result);
                        //System.out.print(" ");
                    }
                }
                j++;
                //System.out.println("a" + j + "\n");
            }
        }catch(Exception e){
            e.printStackTrace();;
        }
        return M1;
    }
    
    public static double[][] ReadMatrixCsv2(String name, int row, int col){
        double [][]M1 = new double[row][col];
        //System.out.println(name);
        try{
            BufferedReader reader = new BufferedReader(new FileReader(name));
            String line = null;
            int j = 0;
            while((line = reader.readLine()) != null){
               	M1[j][col - 1] = Double.parseDouble(line + ".0");
                //System.out.println(line);
                j++;
                //System.out.println();
            }
            //System.out.println("line number" + j);
        }catch(Exception e){
            e.printStackTrace();;
        }
        return M1;
    }
    
    public static void CopyMatrix(double[][] a, double[][] b){
        for(int i = 0; i < b.length; i++){
            for(int j = 0; j < b[i].length; j++){
                a[i][j] = b[i][j];
            }
        }
    }

	public static int findLength(Matrix t, int iter){
		int result = 0;
		double [][]temp = t.getArray();
		for(int i = 0; i < 242; i++){
			if(temp[iter][i] != 0)
				result++;
				//System.out.print(temp[iter][i] + " ");
		}
		return result;
	}
	
	public static double find_abs_max(double[][] m){
		for(int i = 0; i < m.length; i++){
			for(int j = 0; j < m[0].length; j++){
				m[i][j] = Math.abs(m[i][j]);
			}
		}
		double largest = Double.NEGATIVE_INFINITY;
		for(int i = 0; i < m.length; i++){
			for(int j = 0; j < m[0].length; j++){
				if(m[i][j] > largest)
					largest = m[i][j];
			}
		}
		return largest;
	}

	public static void Pause(){
		System.out.println("Press any key to continue");
		try{
			System.in.read();
		}catch(Exception e){

		}
	}

    public static void LoadData(){
        //read string name csv
        System.out.println("read the proteinName file");
        ReadNameFile("./data/proteinName.csv", proteinName);
		//PrintNameList(proteinName);
		System.out.println("read the drugName file");
        ReadNameFile("./data/drugName.csv", drugName);
        
        //read Rp matrix
        System.out.println("read the Rp file");
        Rp = new double[242][1];
        CopyMatrix(Rp, ReadMatrixCsv("./data/Rp.csv", 242, 1));
        Rp_M = new Matrix(Rp);

        //read Expression matrix
        System.out.println("read the Expression file");
        Expression = new double[242][1];
        CopyMatrix(Expression, ReadMatrixCsv2("./data/Exp.csv", 242, 1));
       	Expression_M = new Matrix(Expression);

        //read adj matrix
        System.out.println("read the adj file");
        adj = new double[242][242];
        CopyMatrix(adj, ReadMatrixCsv("./data/adj.csv", 242, 242));
        adj_M = new Matrix(adj);
		//PrintMatrix(adj_M);
		//System.out.println(adj_M.getArray().length + " " + adj_M.getArray()[0].length);

		//read drugVector matrix
        System.out.println("read the drugVector file");
        drugVector = new double[242][97];
        CopyMatrix(drugVector, ReadMatrixCsv("./data/drugVector.csv", 242, 97));
    	drugVector_M = new Matrix(drugVector);
	}
   
	public static void computePETOneDrug(Matrix drugVector_i){
		Matrix proteinWeight_M = new Matrix(242, 1);
		//proteinWeight_M = Rp_M;
		CopyMatrix(Rp, ReadMatrixCsv("./data/Rp.csv", 242, 1));
        Rp_M = new Matrix(Rp);
		double [][]Rp1 = Rp_M.getArray();
		for(int x = 0; x < 242; x++){
			for(int y = 0; y < 1; y++){
				proteinWeight_M.set(x, y, Rp1[x][y]);
			}
		}

		CopyMatrix(adj, ReadMatrixCsv("./data/adj.csv", 242, 242));
		adj_M = new Matrix(adj);
		double [][]absAdj1 = new double[242][242];
		for(int i = 0; i < 242; i++){
			for(int j = 0; j < 242; j++){
				absAdj1[i][j] = Math.abs(adj_M.get(i, j));
			}
		}
		//Matrix absAdj_M = AbsMatrix(adj_M);
		Matrix absAdj_M = new Matrix(absAdj1);

		/*double bb = 0;
		for(int x = 0; x < 10; x++){
			System.out.print(x + " line ");
			for(int y = 0; y < 242; y++){
				if(adj_M.get(x, y) != 0){
					bb++;
					System.out.print(adj_M.get(x, y) + " ");
				}
			}
			System.out.println();
		}
		System.out.println("test1" + bb);
		Pause();*/

		
		for(int i = 0; i < proteinName.size(); i++){
			int sizeDownStream = findLength(adj_M, i);
			if(sizeDownStream >= 1)
				Rp_M.set(i, 0, (double)(sizeDownStream) );
			else{
				Rp_M.set(i, 0, 1);
				sizeDownStream = 1;
			}
			if(sizeDownStream > 0){
				for(int k = 0; k < 242; k++){
					adj_M.set(i, k, adj_M.get(i, k) / (double)(sizeDownStream) );	
				}
			}
			//System.out.println(sizeDownStream);
		}	

		/*double bb = 0;
		for(int x = 0; x < 242; x++){
			System.out.print(x + " line ");
			for(int y = 0; y < 1; y++){
				if(proteinWeight_M.get(x, y) != 0){
					bb++;
					System.out.print(proteinWeight_M.get(x, y) + " ");
				}
			}
			System.out.println();
		}
		System.out.println(bb);
		Pause();*/
		//intialization
		//numOfPro = proteinName.size();
		//~ if the element is zero, it become 1
		//if the element is not zero, it become 0
		//initialize the target matrix to drugVector using ~
		target = new double[242][1];
		for(int i = 0; i < 242; i++){
			if(drugVector_i.get(i, 0) > 0 || drugVector_i.get(i, 0) < 0){
				target[i][0] = 1;
				//System.out.println(i);
			}else{
				target[i][0] = 0;
			}
		}
		Matrix target_M = new Matrix(target);
		//System.out.println(drugVector_M.get(116, 0));
		//now, previous layer and current layer are represented by binary vector
		Matrix previousLayer_M = target_M;

		Matrix absAdj_M_transpose = absAdj_M.transpose();

		Matrix currentLayer_M = CreateSignMatrix(absAdj_M_transpose.times(previousLayer_M));

		/*Matrix ccc = absAdj_M_transpose;
		double aa = 0;
		for(int x = 0; x < 242; x++){
			System.out.print(x + " line ");
			for(int y = 0; y < 1; y++){
				if(currentLayer_M.get(x, y) != 0){
					aa++;
					System.out.print(currentLayer_M.get(x, y) + " ");
				}
			}
			System.out.println();
		}
		System.out.println(aa);
		Pause();*/

		//initialize s score, which is drug-protein interaction
		scoreArray_M = new Matrix(numOfPro, maxIteration, 0.00);
		for(int i = 0; i < numOfPro; i++){
			scoreArray_M.set(i, 0, drugVector_i.get(i, 0) * Rp[i][0]);
			//if(scoreArray_M.get(i, 0) > 0)
			//	System.out.println(scoreArray_M.get(i, 0) + " " + i);
		}
		//System.out.println(scoreArray_M.getRowDimension() + " " + scoreArray_M.getColumnDimension());
		
		double [][]currentlayer = new double[numOfPro][1];
		currentlayer = currentLayer_M.getArray();

		/*double aa = 0, bb = 0;
		for(int x = 0; x < 242; x++){
			System.out.print(x + " line ");
			for(int y = 0; y < 1; y++){
				if(currentLayer_M.get(x, y) != 0){
					aa++;
					System.out.print(currentLayer_M.get(x, y) + " ");
				}
				if(currentlayer[x][y] != 0){
					bb++;
				}
			}
			System.out.println();
		}
		System.out.println(aa + " " + bb);
		Pause();*/



		Matrix previousLayer_M1 = previousLayer_M;
		//399 iteration
		for(int i = 1; i < maxIteration; i++){
			/*double aa = 0;
			for(int x = 0; x < 11; x++){
				System.out.print(x + " line ");
				for(int y = 0; y < 1; y++){
					if(previousLayer_M1.get(x, y) != 0){
						aa++;
						System.out.print(previousLayer_M1.get(x, y) + " ");
					}
				}
				System.out.println();
			}
			System.out.println(aa);
			Pause();*/
			//move the column of the scoreArray backward
			for(int j = 0; j < numOfPro; j++)
				scoreArray_M.set(j, i, scoreArray_M.get(j, i - 1));

			//the first term:a vector, entry (1-d)*c for component in currentlayer, otherwise 0
			Matrix leftTermVector_M = new Matrix(numOfPro, 1, 0);
			for(int x = 0; x < numOfPro; x++){
				leftTermVector_M.set(x, 0, (1 - sigma) * scoreArray_M.get(x, 0) * currentLayer_M.get(x, 0));
				//if(x == 115)
				//	System.out.println((1 - sigma) * scoreArray_M.get(x, 0) + " " + currentLayer_M.get(x, 0));
			}

			/*double aa = 0;
			for(int x = 0; x < 242; x++){
				System.out.print(x + " line ");
				for(int y = 0; y < 1; y++){
					if(leftTermVector_M.get(x, y) != 0){
						aa++;
						System.out.print(leftTermVector_M.get(x, y) + " ");
					}
				}
				System.out.println();
			}
			System.out.println(aa);
			Pause(); */

			//the right term, nodes not belonging to the currentLayer is not updated, therefore create the identity matrix beforehand
			Matrix adjThisTurn_M = CreateEyeMatrix(numOfPro);
			
			//System.out.println(adjThisTurn_M.get(0, 41) + " " + adjThisTurn_M.get(0, 43) );
			//now, for the currentLayer row in adjThisTurn, change the row according to adj matrix
			/*double aa = 0;
			for(int x = 0; x < 242; x++){
				for(int y = 0; y < 1; y++){
					if(currentLayer_M.get(x, y) != 0){
						aa++;
						//System.out.print(currentLayer_M.get(x, y) + " line" + x + " " );
					}
				}
				//System.out.println();
			}
			System.out.println(aa);
			Pause();*/
			for(int j = 0; j < numOfPro; j++){
				/*double bb = 0;
				for(int x = 0; x < 242; x++){
					//System.out.print(x + " line ");
					for(int y = 0; y < 1; y++){
						if(currentLayer_M.get(x, y) != 0){
							bb++;
							//System.out.print(currentLayer_M.get(x, y) + " ");
						}
					}
					//System.out.println();
				}
				System.out.println(bb);*/
				//Pause();

				//System.out.println(adjThisTurn_M.get(0, 41) + " " + adjThisTurn_M.get(0, 43) );
				Matrix previousLayer_M_transpose = previousLayer_M1.transpose();
				Matrix b_vec_M = new Matrix(1, 242, 0);
				for(int k = 0; k < 242; k++){
					b_vec_M.set(0, k, Math.pow(boostingFactor, previousLayer_M_transpose.get(0, k)));

				}
					
				/*if(i == 3 && j == 22) //value should be 1.5
					if(b_vec_M.get(0, 6) == 1.5)
						System.out.println("b_vec_M value right");
					else
						System.out.println("b_vec_M value wrong " + b_vec_M.get(0, 6));*/
				//get adj(:, i)
				//System.out.println("read the adj file");
		        //adj = new double[242][242];
		        //CopyMatrix(adj, ReadMatrixCsv("./data/adj.csv", 242, 242));
		        adj_M = new Matrix(adj);
				
				Matrix adj_col_M = adj_M.getMatrix(0, 241,  j, j).transpose();
				//if(j == 22)
				//	System.out.println(adj_col_M.get(0, 6) );

				/*for(int x = 0; x < 1; x++){
					for(int y = 6; y < 7; y++){
						if(i == 3 && j == 7)
							System.out.println(adj_col_M.get(x, y));
					}

				}*/

				/*double aa = 0;
				for(int x = 0; x < 1; x++){
					System.out.print(x + " line ");
					for(int y = 0; y < 7; y++){
						//if(vvv.get(x, y) != 0){
						//	aa++;
							System.out.print(adj_col_M.get(x, y) + " ");
						//}
					}
					System.out.println();
				}
				//System.out.println(aa);
				Pause();*/
				double [][]adj_col = new double[1][242];
				adj_col = adj_col_M.getArray();
				Matrix adj_col_MMM = new Matrix(adj_col);
				Matrix newRow_M = b_vec_M.arrayTimes(adj_col_MMM.times(sigma));
				//Matrix newRow_M = b_vec_M.arrayTimes(adj_col_M.times(sigma));
				/*if(i == 3 && j == 20){
					for(int x = 0; x < 1; x++){
						//System.out.print(x + " line ");
						for(int y = 0; y < 10; y++){
							System.out.print(newRow_M.get(x,y) + " ");
						}
						//System.out.println();
					}
				}*/
				
				if(currentLayer_M.get(j, 0) == 1){
					//adjThisTurn_M.setMatrix(j, j, 0, 241, newRow_M.times(currentLayer_M.get(j, 0)) );
					Matrix xxxx = newRow_M.times(currentLayer_M.get(j, 0));
					for(int w = 0; w < 242; w++){
						adjThisTurn_M.set(j, w, xxxx.get(0, w));
						//System.out.println(adjThisTurn_M.get(j, w) + " " + xxxx.get(0, w));
						//ERROR
						//adjthisturn is not return the right result, seems the set function is not right
					}
					/*double aa = 0;
					for(int x = 0; x < 1; x++){
						for(int y = 0; y < 242; y++){
							if(xxxx.get(x, y) != 0){
								aa++;
								System.out.print(xxxx.get(x, y) + " line" + y + " " );
							}
						}
						System.out.println();
					}
					System.out.println(aa);
					Pause();*/
				}
				/*if(i == 2 &&  j == 241)
					System.out.print(adjThisTurn_M.get(12, 6) + " ");*/
				/*double aa = 0;
				for(int x = 0; x < 242; x++){
					//System.out.print(x + " line ");
					for(int y = 0; y < 242; y++){
						if(adjThisTurn_M.get(x, y) != 0){
							aa++;
							//System.out.print(adjThisTurn_M.get(x, y) + " " + x + " line \n");
						}
					}
					//System.out.println();
				}
				System.out.println(aa);
				Pause();*/
				//System.out.println(currentLayer_M.get(58, 0) +  " " + currentLayer_M.get(59, 0) + " " + currentLayer_M.get(60, 0));
				//Pause();
				/*/double aa = 0;
				for(int x = 0; x < 242; x++){
					System.out.print(x + " line ");
					for(int y = 0; y < 1; y++){
						if(currentLayer_M.get(x, y) != 0){
							aa++;
							System.out.print(currentLayer_M.get(x, y) + " ");
						}
					}
					System.out.println();
				}
				System.out.println(aa);
				Pause();*/
				
			}
			/*double aa = 0;
			for(int x = 0; x < 15; x++){
				System.out.print(x + " line ");
				for(int y = 0; y < 15; y++){
					if(adjThisTurn_M.get(x, y) != 0){
						aa++;
						System.out.print(adjThisTurn_M.get(x, y) + " ");
					}
				}
				System.out.println();
			}
			System.out.println(aa);
			Pause();*/
			
			//scoreArray
			/*double aa = 0;
			for(int x = 0; x < 242; x++){
				System.out.print(x + " line ");
				for(int y = 0; y < 150; y++){
					if(scoreArray_M.get(x, y) != 0){
						aa++;
						System.out.print(scoreArray_M.get(x, y) + " ");
					}
				}
				System.out.println();
			}
			System.out.println(aa);
			Pause();*/

			//ERROR the 53 line of the scoreArray in the 3rd iteration is wrong
			
			//leftterm is wrong
			/*double aa = 0;
			for(int x = 0; x < 242; x++){
				System.out.print(x + " line ");
				for(int y = 0; y < 1; y++){
					if(leftTermVector_M.get(x, y) != 0){
						aa++;
						System.out.print(leftTermVector_M.get(x, y) + " ");
					}
				}
				System.out.println();
			}
			System.out.println(aa);
			Pause(); */

			Matrix temp_r1_M = leftTermVector_M.plus(adjThisTurn_M.times(scoreArray_M.getMatrix(0, 241, i - 1, i - 1)));
			//Matrix temp1 = scoreArray_M.getMatrix(0, 241, i - 1, i - 1);
			scoreArray_M.setMatrix(0, 241, i, i, temp_r1_M);
			
			//ERROR from i > 4
			/*double aa = 0;
			for(int x = 0; x < 242; x++){
				System.out.print(x + " line ");
				for(int y = 0; y < 242; y++){
					if(adjThisTurn_M.get(x, y) != 0){
						aa++;
						System.out.print(adjThisTurn_M.get(x, y) + " ");
					}
				}
				System.out.println();
			}
			System.out.println(aa);
			Pause();*/

			//System.out.println(adjThisTurn_M.get(121, 0) );
			/*for(int y = 0; y < 242; y++)
				if(leftTermVector_M.get(y, 0) > 0)
					System.out.println(leftTermVector_M.get(y, 0));*/
			
			//System.out.println(scoreArray_M.get(6, 1) + " " +  scoreArray_M.get(6, 2) + " " +  scoreArray_M.get(6, 3));
			//update layer
			previousLayer_M1 = currentLayer_M;
			//System.out.println(previousLayer_M.get(6, 0) + " " + previousLayer_M.get(52, 0) );
			//The currentLayer here is right but when enter the next interation, it is init wrong and go back to the init state
			/*double aa = 0;
			for(int x = 0; x < 242; x++){
				//System.out.print(x + " line ");
				for(int y = 0; y < 1; y++){
					if(currentLayer_M.get(x, y) != 0){
						aa++;
						//System.out.print(currentLayer_M.get(x, y) + " ");
					}
				}
				//System.out.println();
			}
			System.out.println(aa);
			Pause();*/
			currentLayer_M = CreateSignMatrix(absAdj_M_transpose.times(previousLayer_M1));
			/*double aa = 0;
			for(int x = 0; x < 242; x++){
				//System.out.print(x + " line ");
				for(int y = 0; y < 1; y++){
					if(currentLayer_M.get(x, y) != 0){
						aa++;
						//System.out.print(currentLayer_M.get(x, y) + " ");
					}
				}
				//System.out.println();
			}
			System.out.println(aa);
			Pause();*/
		}

		/*double aa = 0;
		for(int x = 0; x < 10; x++){
			System.out.print(x + " line ");
			for(int y = 0; y < 10; y++){
				if(scoreArray_M.get(x, y) != 0){
					aa++;
					System.out.print(scoreArray_M.get(x, y) + " ");
				}
			}
			System.out.println();
		}
		System.out.println(aa);
		Pause();*/

		/*double aa = 0;
		for(int x = 0; x < 242; x++){
			System.out.print(x + " line ");
			for(int y = 0; y < 1; y++){
				if(currentLayer_M.get(x, y) != 0){
					aa++;
					System.out.print(currentLayer_M.get(x, y) + " ");
				}
			}
			System.out.println();
		}
		System.out.println(aa);
		Pause();*/
		//PrintMatrix(scoreArray_M.getMatrix(0, 4, 0, 4));
		Matrix sMat_M = scoreArray_M.getMatrix(0, 241, 149, 149);
		
		/*for(int i = 0; i < 10; i++){
			for(int j = 0; j < 1; j++){
				System.out.print(sMat_M.get(i, j) + "   ");
			}
			System.out.println();
		}*/
		for(int i = 0; i < 242; i++){
			Matrix temp_r2 = scoreArray_M.getMatrix(i, i, 0, 149);
			double maxSStream = find_abs_max(temp_r2.getArray());
			//System.out.println(maxSStream);
			if( (Math.abs(scoreArray_M.get(i, 149)) / maxSStream) < 0.1){
				sMat_M.set(i, 0, 0);
			}
		}
		/*for(int i = 0; i < 10; i++){
			for(int j = 0; j < 1; j++){
				System.out.print(sMat_M.get(i, j) + "   ");
			}
			System.out.println();
		}*/

		//compute drug ranking score after the given iteration
		double nume = 0;
		double deno = 0;
		for(int i = 0; i < proteinName.size(); i++){
			double val = 0.0;
			//System.out.println(sMat_M.get(i, 0) + " " + Expression[i][0]);
			if(Expression[i][0] != 0){
				if(sMat_M.get(i, 0) * Expression[i][0] > 0)
					val = 1;
				else if(sMat_M.get(i, 0) * Expression[i][0] == 0)
					val = 0;
				else if(sMat_M.get(i, 0) * Expression[i][0] < 0)
					val = -1;
				nume -= proteinWeight_M.get(i, 0) * val;
				deno += proteinWeight_M.get(i, 0);
			}
		}

		if(deno != 0){
			PETScoreMatrix = nume / deno;	
		}else{
			PETScoreMatrix = 0;
		}
		System.out.println("The score is " + PETScoreMatrix);
	}
	
    public static void main(String argv[]){
        System.out.println("PETS alroithm begin\n");
        
		//load the matrix in the data folder
        LoadData();
        
		//algorithm begin
		//get the value of  the drugPro martix
		//System.out.print(drugName.size() + " " + proteinName.size());
		double [][]drugpro = new double[drugName.size()][proteinName.size()];
		InitDoubleArray(drugpro);
		Matrix drugPro_M = new Matrix(drugpro);
	
		//reference the computePETOneDrug function
		int useProteinDegree = 0;
		ArrayList<Double> drugRankingScore = new ArrayList<Double>();
		for(int i = 0; i < drugName.size(); i++){
			//due to the useProteinDegree is never used, I delete it
			//computer PET one drug
			Matrix drugVector_i;
			drugVector_i = drugVector_M.getMatrix(0, 241,  i, i);
			computePETOneDrug(drugVector_i);
			drugRankingScore.add(i, PETScoreMatrix);	
			//Pause();
		}
		
		System.out.println("The score is " + PETScoreMatrix);
		
    }
}

























