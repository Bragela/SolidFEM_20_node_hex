using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using LA = MathNet.Numerics.LinearAlgebra;
using CSparse;
using CSD = CSparse.Double;
using CSparse.Double.Factorization;
using CSparse.Storage;

using Rhino.Geometry;

namespace SolidFEM.Classes
{
    /// <summary>
    /// A class for all the matrix methods needed in the finite element analysis. 
    /// </summary>
    public static class FEM_Matrices
    {
        /// <summary>
        /// Construct global stiffness matrix by assembling element stiffness matrices.
        /// </summary>
        /// <returns> Global stiffness matrix. </returns>
        public static double[,] GlobalStiffnessCSparse(ref List<Element> elements, int numNode, Material material, ref FEMLogger logger)
        {
            Stopwatch timer = new Stopwatch();

            //Create empty K matrix
            double[,] kArray = new double[numNode * 3, numNode * 3];
                
            //Loop through all elements
            foreach (Element element in elements)
            {
                List<int> con = element.Connectivity; // get the connectivity of each element

                //Iterate over the connectivity indices
                var kAndB = CalculateElementMatrices(element, material, ref logger, "Full");
                LA.Matrix<double> K_local = kAndB.Item1;
                element.LocalB = kAndB.Item2;
                
                //Loop nodes of elements
                for (int i = 0; i < con.Count; i++)
                {
                    for (int j = 0; j < con.Count; j++)
                    {
                        //Loop relevant local stiffness contribution. Add local stiffness matrix to correct index in global K
                        for (int dofRow = 0; dofRow < 3; dofRow++)
                        {
                            for (int dofCol = 0; dofCol < 3; dofCol++)
                            {
                                int lID1 = 3 * i + dofRow;
                                int lID2 = 3 * j + dofCol;
                                int gID1 = 3 * con[i] + dofRow;
                                int gID2 = 3 * con[j] + dofCol;

                                kArray[gID1, gID2] += K_local[lID1, lID2];
                            }
                        }
                    }
                }
            }

            return kArray;
        }

        /// <summary>
        /// Calculate element stifness matrix and element strain matrix.
        /// </summary>
        /// <returns> Element stiffness and strain matrix.</returns>
        public static Tuple<LA.Matrix<double>, List<LA.Matrix<double>>> CalculateElementMatrices(Element element, Material material, ref FEMLogger logger, String intType)
        {
            //Summary: calculate local K and B matrix
            int roundrecisionBMatrix = 6;
            int rpb = roundrecisionBMatrix;
            
            //Material
            LA.Matrix<double> C = material.GetMaterialConstant();

            //Create empty local stiffness matrix
            int numElementNodes = element.Nodes.Count;
            LA.Matrix<double> K_local = LA.Matrix<double>.Build.Dense(3 * numElementNodes, 3 * numElementNodes);

            //Create empty local deformation matrix
            List<LA.Matrix<double>> B_local = new List<LA.Matrix<double>>();

            //Global coordinates of the nodes of the actual element
            LA.Matrix<double> globalCoordinates = LA.Matrix<double>.Build.Dense(numElementNodes, 3);
            List<Point3d> localCoordinates = FEM_Utility.LocalCartesianCoordinates(element);

            for (int i = 0; i < numElementNodes; i++)
            {
                globalCoordinates[i, 0] = Math.Round(element.Nodes[i].Coordinate.X, rpb); // column of x coordinates
                globalCoordinates[i, 1] = Math.Round(element.Nodes[i].Coordinate.Y, rpb); // column of y coordinates
                globalCoordinates[i, 2] = Math.Round(element.Nodes[i].Coordinate.Z, rpb); // column of z coordinates
            }

            // Different methods for Hex8, Hex20 and Tet4. Tet4 doesn't need gauss integration because B and Jacobian are constant!

            if (element.Type == "Hex8")
            {
                //Numerical integration
                //LA.Matrix<double> gaussNodes = FEM.GetNaturalCoordinate((double)Math.Sqrt((double)1 / (double)3), 3);
                var gaussCoordinates = FEM_Utility.GetGaussPointMatrix(2, element.Type); // by defaul we have a 2x2x2 integration of Hex8 element
                List<double> pointJacobians = new List<double>(); // list to evaluate the pointwise jacobians
                for (int n = 0; n < gaussCoordinates.RowCount; n++)  // loop gauss nodes
                {
                    // Substitute the natural coordinates into the symbolic expression
                    var r = gaussCoordinates.Row(n)[0];
                    var s = gaussCoordinates.Row(n)[1];
                    var t = gaussCoordinates.Row(n)[2];

                    // Partial derivatives of the shape functions
                    //LA.Matrix<double> shapeFunctionsDerivatedNatural = FEM.DerivateWithNatrualCoordinates(r, s, t, 3);
                    var partialDerivatives = FEM_Utility.PartialDerivateShapeFunctions(r, s, t, element.Type);

                    // Calculate Jacobian matrix
                    var jacobianMatrix = partialDerivatives.Multiply(globalCoordinates);

                    LA.Matrix<double> shapeFuncDerivatedCartesian = (jacobianMatrix.Inverse()).Multiply(partialDerivatives);

                    // calculate B with CSparse




                    double jacobianDeterminant = jacobianMatrix.Determinant(); // To do: Find a way to evaluate this...
                    pointJacobians.Add(jacobianDeterminant);

                    if (jacobianDeterminant < 0) { logger.AddWarning("Negativ jacobian determeninant"); }
                    int dimRowB = 6;

                    // establish the B-matrix
                    LA.Matrix<double> B_i = LA.Double.DenseMatrix.Build.Dense(dimRowB, 3 * numElementNodes);

                    for (int i = 0; i < numElementNodes; i++)
                    {

                        // with the shape functions derivated with respect to the cartesian coordinates the rotated and unrotated element vectors are not the same... This is the correct one according to the formulas
                        var B_i_sub = LA.Double.DenseMatrix.Build.DenseOfRowMajor(6, 3, new double[] {
                            shapeFuncDerivatedCartesian.Row(0)[i], 0, 0,
                            0, shapeFuncDerivatedCartesian.Row(1)[i], 0,
                            0, 0, shapeFuncDerivatedCartesian.Row(2)[i],
                            shapeFuncDerivatedCartesian.Row(1)[i], shapeFuncDerivatedCartesian.Row(0)[i], 0,
                            shapeFuncDerivatedCartesian.Row(2)[i], 0, shapeFuncDerivatedCartesian.Row(0)[i],
                            0, shapeFuncDerivatedCartesian.Row(2)[i], shapeFuncDerivatedCartesian.Row(1)[i]
                            });

                        B_i.SetSubMatrix(0, i * 3, B_i_sub);

                    }

                    B_local.Add(B_i);
                    var k_i = (B_i.Transpose()).Multiply(C.Multiply(B_i)).Multiply(jacobianDeterminant);

                    K_local.Add(k_i, K_local);
                }
            }
            else if (element.Type == "Tet4")
            {
                var partialDerivatives = FEM_Utility.PartialDerivateShapeFunctions(1, 1, 1, "Tet4");    //Random coordinates (1,1,1) because the method requires coordinate inputs

                // Calculate Jacobian matrix
                var jacobianMatrix = partialDerivatives.Multiply(globalCoordinates);

                // Calculate B - LA.Matrix
                LA.Matrix<double> shapeFuncDerivatedCartesian = (jacobianMatrix.Inverse()).Multiply(partialDerivatives);

                int dimRowB = 6;
                LA.Matrix<double> B_i = LA.Double.DenseMatrix.Build.Dense(dimRowB, 3 * numElementNodes);

                for (int i = 0; i < numElementNodes; i++)
                {

                    // with the shape functions derivated with respect to the cartesian coordinates the rotated and unrotated element vectors are not the same... This is the correct one according to the formulas
                    var B_i_sub = LA.Double.DenseMatrix.Build.DenseOfRowMajor(6, 3, new double[] {
                            shapeFuncDerivatedCartesian.Row(0)[i], 0, 0,
                            0, shapeFuncDerivatedCartesian.Row(1)[i], 0,
                            0, 0, shapeFuncDerivatedCartesian.Row(2)[i],
                            shapeFuncDerivatedCartesian.Row(1)[i], shapeFuncDerivatedCartesian.Row(0)[i], 0,
                            shapeFuncDerivatedCartesian.Row(2)[i], 0, shapeFuncDerivatedCartesian.Row(0)[i],
                            0, shapeFuncDerivatedCartesian.Row(2)[i], shapeFuncDerivatedCartesian.Row(1)[i]
                            });


                    B_i.SetSubMatrix(0, i * 3, B_i_sub);

                }
                B_local.Add(B_i);
                // Get volume of Tetrahedra
                VolumeMassProperties vmp = VolumeMassProperties.Compute(element.ElementMesh);
                double V = vmp.Volume;

                var k_i = V * (B_i.Transpose()).Multiply(C.Multiply(B_i));

                K_local.Add(k_i, K_local);
            }
            else if (element.Type == "Hex20")
            {
                if (intType == "Full")
                {
                    //Numerical integration

                    int order = 3;
                    var gaussCoordinates = FEM_Utility.GetGaussPointMatrix(order, element.Type);    // by default we have a 3x3x3 integration of Hex20 element
                    List<double> pointJacobians = new List<double>();                               // list to evaluate the pointwise jacobians
                    for (int n = 0; n < gaussCoordinates.RowCount; n++)                             // loop gauss nodes
                    {
                        // Substitute the natural coordinates into the symbolic expression
                        var r = gaussCoordinates.Row(n)[0];
                        var s = gaussCoordinates.Row(n)[1];
                        var t = gaussCoordinates.Row(n)[2];

                        // Partial derivatives of the shape functions
                        var partialDerivatives = FEM_Utility.PartialDerivateShapeFunctions(r, s, t, element.Type);

                        // Calculate Jacobian matrix
                        var jacobianMatrix = partialDerivatives.Multiply(globalCoordinates);

                        // Derivated shape function in cartesian coordinates
                        LA.Matrix<double> shapeFuncDerivatedCartesian = (jacobianMatrix.Inverse()).Multiply(partialDerivatives);

                        double jacobianDeterminant = jacobianMatrix.Determinant();
                        pointJacobians.Add(jacobianDeterminant);

                        if (jacobianDeterminant < 0) 
                        { logger.AddWarning("Negativ jacobian determeninant"); }
                        int dimRowB = 6;


                        //Establish the B-matrix
                        LA.Matrix<double> B_i = LA.Double.DenseMatrix.Build.Dense(dimRowB, 3 * numElementNodes);

                        for (int i = 0; i < numElementNodes; i++)
                        {
                            var B_i_sub = LA.Double.DenseMatrix.Build.DenseOfRowMajor(6, 3, new double[] {
                            shapeFuncDerivatedCartesian.Row(0)[i], 0, 0,
                            0, shapeFuncDerivatedCartesian.Row(1)[i], 0,
                            0, 0, shapeFuncDerivatedCartesian.Row(2)[i],
                            shapeFuncDerivatedCartesian.Row(1)[i], shapeFuncDerivatedCartesian.Row(0)[i], 0,
                            shapeFuncDerivatedCartesian.Row(2)[i], 0, shapeFuncDerivatedCartesian.Row(0)[i],
                            0, shapeFuncDerivatedCartesian.Row(2)[i], shapeFuncDerivatedCartesian.Row(1)[i]
                            });
                            B_i.SetSubMatrix(0, i * 3, B_i_sub);
                        }

                        B_local.Add(B_i);

                        //Weights for numerical integration

                        double w_r = 0.555556; double w_s = 0.5555556; double w_t = 0.555556;

                        if (r == 0)
                        {
                            w_r = 0.888889;
                        }
                        if (s == 0)
                        {
                            w_s = 0.888889;
                        }
                        if (t == 0)
                        {
                            w_t = 0.888889;
                        }

                        //Calculate local K matrix with gaussian numerical integration
                        var k_i = (B_i.Transpose()).Multiply(C.Multiply(B_i)).Multiply(jacobianDeterminant) * w_r * w_s * w_t;

                        K_local.Add(k_i, K_local);
                    }
                    
                }
                else if (intType == "Reduced")
                {
                    //Numerical integration

                    int order = 2;
                    var gaussCoordinates = FEM_Utility.GetGaussPointMatrix(order, element.Type);    // with reduced integration we have a 2x2x2 integration of Hex20 element 
                    List<double> pointJacobians = new List<double>();                               // list to evaluate the pointwise jacobians
                    for (int n = 0; n < gaussCoordinates.RowCount; n++)                             // loop gauss nodes
                    {
                        // Substitute the natural coordinates into the symbolic expression
                        var r = gaussCoordinates.Row(n)[0];
                        var s = gaussCoordinates.Row(n)[1];
                        var t = gaussCoordinates.Row(n)[2];

                        // Partial derivatives of the shape functions
                        var partialDerivatives = FEM_Utility.PartialDerivateShapeFunctions(r, s, t, element.Type);

                        // Calculate Jacobian matrix
                        var jacobianMatrix = partialDerivatives.Multiply(globalCoordinates);

                        // Calculate B - LA.Matrix

                        LA.Matrix<double> shapeFuncDerivatedCartesian = (jacobianMatrix.Inverse()).Multiply(partialDerivatives);

                        // calculate B with CSparse

                        double jacobianDeterminant = jacobianMatrix.Determinant(); // To do: Find a way to evaluate this...
                        pointJacobians.Add(jacobianDeterminant);

                        if (jacobianDeterminant < 0) { logger.AddWarning("Negativ jacobian determeninant"); }
                        int dimRowB = 6;


                        // establish the B-matrix
                        LA.Matrix<double> B_i = LA.Double.DenseMatrix.Build.Dense(dimRowB, 3 * numElementNodes);

                        for (int i = 0; i < numElementNodes; i++)
                        {
                            var B_i_sub = LA.Double.DenseMatrix.Build.DenseOfRowMajor(6, 3, new double[] {
                            shapeFuncDerivatedCartesian.Row(0)[i], 0, 0,
                            0, shapeFuncDerivatedCartesian.Row(1)[i], 0,
                            0, 0, shapeFuncDerivatedCartesian.Row(2)[i],
                            shapeFuncDerivatedCartesian.Row(1)[i], shapeFuncDerivatedCartesian.Row(0)[i], 0,
                            shapeFuncDerivatedCartesian.Row(2)[i], 0, shapeFuncDerivatedCartesian.Row(0)[i],
                            0, shapeFuncDerivatedCartesian.Row(2)[i], shapeFuncDerivatedCartesian.Row(1)[i]
                            });
                            B_i.SetSubMatrix(0, i * 3, B_i_sub);
                        }

                        B_local.Add(B_i);

                        //Weights for numerical integration

                        double w_r = 1; double w_s = 1; double w_t = 1;

                        //Calculate local K matrix with gaussian numerical integration
                        var k_i = (B_i.Transpose()).Multiply(C.Multiply(B_i)).Multiply(jacobianDeterminant) * w_r * w_s * w_t;

                        K_local.Add(k_i, K_local);
                    }
                }
            }
            return Tuple.Create(K_local, B_local);
        }

        public static LA.Matrix<double> CalculateGlobalStiffnessMatrix(List<Element> elements, int numNode, Material material, ref FEMLogger logger)
        {

            // create stiffness matrix
            LA.Matrix<double> K_global = LA.Matrix<double>.Build.Dense(numNode * 3, numNode * 3);
            foreach (Element element in elements)
            {

                List<int> con = element.Connectivity;


                LA.Matrix<double> K_local = FEM_Matrices.CalculateElementMatrices(element, material, ref logger, "Full").Item1;

                // loop nodes of elements
                for (int i = 0; i < con.Count; i++)
                {
                    for (int j = 0; j < con.Count; j++)
                    {
                        // loop relevant local stiffness contribution
                        for (int dofRow = 0; dofRow < 3; dofRow++)
                        {
                            for (int dofCol = 0; dofCol < 3; dofCol++)
                            {
                                K_global[3 * con[i] + dofRow, 3 * con[j] + dofCol] += +K_local[3 * i + dofRow, 3 * j + dofCol];
                            }
                        }
                    }
                }

            }
            return K_global;
        }


        
    }
}
