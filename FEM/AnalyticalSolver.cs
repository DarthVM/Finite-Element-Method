using MathNet.Numerics.LinearAlgebra;
using static System.Math;


namespace FEM;

public class AnalyticalSolver(InitialConditions init)
{
    
    private Vector<double> C(Matrix<double> A, Vector<double> vectorB)
    {
        return A.Solve(vectorB);
    }

    public Func<double, double> Solve()
    {
        
        var disc = init.B * init.B - 4 * init.A * init.C;

        if (disc >= 0)
        {
            var l = new[]
            {
                (-init.B + Math.Sqrt(disc)) / (2 * init.A),
                (-init.B - Math.Sqrt(disc)) / (2 * init.A)
            };

            var (A, vectorB) = DiscriminantPos(l);
            Vector<double> constant = C(A, vectorB);
            return x =>
                constant[0] * Math.Exp(l[0] * x) +
                constant[1] * Math.Exp(l[1] * x) -
                (init.C != 0 ? (double)init.D / init.C : init.D * x / init.B);
        }
        else
        {
            var alpha = -init.B;
            var beta = Math.Sqrt(Math.Abs(disc));


            var (A, vectorB) = DiscriminantNeg(alpha, beta);
            var constant = C(A, vectorB);
            return x =>
                constant[0] * Math.Exp(alpha * x) * Math.Cos(beta * x) +
                constant[1] * Math.Exp(alpha * x) * Math.Sin(beta * x) +
                (double)init.D / (init.C != 0 ? init.C : init.B);

        }
    }

    

    private (Matrix<double>, Vector<double>) DiscriminantPos(double[] lambda)
    {
        Matrix<double> A = Matrix<double>.Build.Dense(2, 2);
        Vector<double> vectorB = Vector<double>.Build.Dense(2);


        switch(init.LeftCondition.Type)
        {

            case ConditionType.First:

                A[0, 0] = Math.Exp(lambda[0] * init.StartX);

                A[0, 1] = Math.Exp(lambda[1] * init.StartX);

                vectorB[0] = init.LeftCondition.Value;

                break;


            case ConditionType.Second:

                A[0, 0] = lambda[0] * Math.Exp(lambda[0] * init.StartX);

                A[0, 1] = lambda[1] * Math.Exp(lambda[1] * init.StartX);

                vectorB[0] = init.LeftCondition.Value;

                break;


            case ConditionType.Third:

                A[0, 0] = Math.Exp(lambda[0] * init.StartX) -
                          lambda[0] * Math.Exp(lambda[0] * init.StartX);


                A[0, 1] = Math.Exp(lambda[1] * init.StartX) -
                          lambda[1] * Math.Exp(lambda[1] * init.StartX);

                vectorB[0] = -(double)init.D / (init.C != 0 ? init.C : init.B) + (double)init.D / (init.C != 0 ? init.C : init.B) * init.StartX;

                break;

        }

        switch (init.RightCondition.Type)
        {

            case ConditionType.First:

                A[1, 0] = Math.Exp(lambda[0] * init.EndX);

                A[1, 1] = Math.Exp(lambda[1] * init.EndX);

                vectorB[1] = init.RightCondition.Value + (double)init.D / (init.C != 0 ? init.C : init.B) * init.EndX;

                break;


            case ConditionType.Second:

                A[1, 0] = lambda[0] * Math.Exp(lambda[0] * init.EndX);

                A[1, 1] = lambda[1] * Math.Exp(lambda[1] * init.EndX);

                vectorB[1] = init.LeftCondition.Value;

                break;


            case ConditionType.Third:

                A[1, 0] = Math.Exp(lambda[0] * init.EndX) -
                          lambda[0] * Exp(lambda[0] * init.EndX);

                A[1, 1] = Exp(lambda[1] * init.EndX) -
                          lambda[1] * Exp(lambda[1] * init.EndX);

                vectorB[1] = 0;

                break;

        }

        return (A, vectorB);
    }




    private (Matrix<double>, Vector<double>) DiscriminantNeg(double alpha, double beta)
    {
        var A = Matrix<double>.Build.Dense(2, 2);
        var vectorB = Vector<double>.Build.Dense(2);


        switch (init.LeftCondition.Type)
        {

            case ConditionType.First:

                A[0, 0] = Math.Exp(alpha * init.StartX) * Math.Cos(beta * init.StartX);

                A[0, 1] = Math.Exp(alpha * init.StartX) * Math.Sin(beta * init.StartX);

                vectorB[0] = init.LeftCondition.Value;

                break;


            case ConditionType.Second:

                A[0, 0] = (alpha + beta) * Math.Exp(alpha * init.StartX) * Math.Cos(beta * init.StartX);

                A[0, 1] = -Math.Exp(alpha * init.StartX) * Math.Sin(beta * init.StartX);

                vectorB[0] = init.LeftCondition.Value;

                break;


            case ConditionType.Third:

                A[0, 0] = (alpha + beta - 1) * Math.Exp(alpha * init.StartX) * Math.Cos(beta * init.StartX);

                A[0, 1] = -2 * Math.Exp(alpha * init.StartX) * Math.Sin(beta * init.StartX);

                vectorB[0] = 0;

                break;

        }

        switch (init.RightCondition.Type)
        {

            case ConditionType.First:

                A[1, 0] = Exp(alpha * init.EndX) * Cos(beta * init.EndX);

                A[1, 1] = Exp(alpha * init.EndX) * Sin(beta * init.EndX);

                vectorB[1] = init.RightCondition.Value;

                break;


            case ConditionType.Second:

                A[1, 0] = (alpha + beta) * Exp(alpha * init.EndX) * Cos(beta * init.EndX);

                A[1, 1] = -Exp(alpha * init.EndX) * Sin(beta * init.EndX);

                vectorB[1] = init.RightCondition.Value;

                break;


            case ConditionType.Third:

                A[1, 0] = (alpha + beta - 1) * Exp(alpha *init.EndX) * Cos(beta * init.EndX);

                A[1, 1] = -2 * Exp(alpha * init.EndX) * Sin(beta *init.EndX);

                vectorB[1] = 0;

                break;

        }

        return (A, vectorB);
    }
}