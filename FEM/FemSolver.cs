using static System.Double;
using System.Diagnostics;
using System.Text;
using Control = MathNet.Numerics.Control;
using MathNet.Numerics.LinearAlgebra;
using MathNet.Numerics.LinearAlgebra.Double;




namespace FEM;

public class FemSolver
{
    private readonly int _a;
    private readonly int _b;
    private readonly int _c;
    private readonly int _d;
    private readonly int _startX;
    private readonly int _endX;
    private readonly EdgeCondition _leftCondition;
    private readonly EdgeCondition _rightCondition;
    
    private int _ElementsNumber;
    private int _PointsNumber;
    private readonly int _length;
    private double[] _grid;
    



    public FemSolver(InitialConditions initialConditions, int elementsNumber
    )
    {
        // Init coefficients
        _a = initialConditions.A;
        _b = initialConditions.B;
        _c = initialConditions.C;
        _d = initialConditions.D;

        _leftCondition = initialConditions.LeftCondition;
        _rightCondition = initialConditions.RightCondition;
        
        
        _startX = initialConditions.StartX;
        _endX = initialConditions.EndX;

        _ElementsNumber = elementsNumber;
        _PointsNumber = elementsNumber + 1;

        _length = (int)Abs(_endX - _startX);
        
        



        Control.UseNativeMKL();
        
    }

    


    private Matrix<double> LinHardMat()
    {


        var lhm = Matrix<double>.Build.Sparse(_PointsNumber, _PointsNumber);
        var last = 0.0;



        for (var i = 0; i < _ElementsNumber; i++)
        {
            Matrix<double> local = Matrix<double>.Build.Dense(2, 2);
            double length = Length(_grid[i], _grid[i + 1]);
            
            if (_a != 0)
            {
                Matrix<double> aMat = DenseMatrix.OfArray(new[,]
                {
                    { 1 / length, - 1 / length },
                    { - 1 / length, 1 / length }
                });
                local = - _a * aMat;
            }
            
            
            
            if (_b != 0)
            {
                Matrix<double> bMat = DenseMatrix.OfArray(new[,]
                {
                    { -0.5, 0.5 },
                    { -0.5, 0.5 }
                });
                local += _b * bMat;
            }

            if (_c != 0)
            {
                Matrix<double> cMat = DenseMatrix.OfArray(new[,]
                {
                    { length / 3.0, length / 6.0 },
                    { length / 6.0, length / 3.0 }
                });
                local += _c * cMat;
            }


            if (i != 0)
            {
                local[0, 0] += last;
                lhm.SetSubMatrix(i, i, local);
            }
            else
                lhm.SetSubMatrix(i, i, local);

            last = local.Diagonal().Last();


        }



        return lhm;



    }


    private Vector<double> LinForceMat()
    {
        Vector<double> lfm = Vector<double>.Build.Dense(_PointsNumber);

        lfm[0] = _d * (Length(_grid[0], _grid[1]) / 2.0);
        lfm[^1] = _d * (Length(_grid[^2], _grid[^1]) / 2.0);

        for (var i = 1; i < _ElementsNumber; i++)
            lfm[i] = _d * (Length(_grid[i - 1], _grid[i]) + Length(_grid[i], _grid[i + 1])) / 2;

        return lfm;
    }




    public Vector<double> LinSolve()
    {
        
        _grid = new double[_PointsNumber];
        _grid[0] = _startX;
        _grid[^1] = _endX;
        
        
        var step = (double)_length / _ElementsNumber;
        for (var i = 1; i < _PointsNumber - 1; i++)
            _grid[i] = _grid[i - 1] + step;


        
        Matrix<double> hlm = LinHardMat();
        Vector<double> flm = LinForceMat();


        


        _leftCondition.Init(ref hlm, ref flm);
        _rightCondition.Init(ref hlm, ref flm);;
        
        
        if (_leftCondition.Type == ConditionType.Third)
        {
            
           return hlm.Svd().Solve(-flm).SubVector(1, _PointsNumber);
           
        }

        return _rightCondition.Type == ConditionType.Third ? 
            hlm.Svd().Solve(-flm).SubVector(0, _ElementsNumber) : 
            hlm.Svd().Solve(-flm);
    }


    public Vector<double> CubSolve()
    {
        _grid = new double[_ElementsNumber * 3 + 1];
        _grid[0] = _startX;
        _grid[^1] = _endX;
        
        
        var step = (double)_length / (_ElementsNumber * 3);
        for (var i = 1; i < _ElementsNumber * 3; i++)
            _grid[i] = _grid[i - 1] + step;
        
        
        
        var chm = CubHardMat();
        var cfm = CubForceMat();

        
        
        
        _leftCondition.Init(ref chm, ref cfm);
        _rightCondition.Init(ref chm, ref cfm);
        
     
        
        if (_leftCondition.Type == ConditionType.Third)
        {
            
            return chm.Svd().Solve(-cfm).SubVector(1, chm.ColumnCount - 1);
        }

        return _rightCondition.Type == ConditionType.Third ? 
            chm.Svd().Solve(-cfm).SubVector(0, chm.ColumnCount - 2) : 
            chm.Svd().Solve(-cfm);
    }


    private Matrix<double> CubHardMat()
    {


        Matrix<double> chm = Matrix<double>.
            Build.Sparse(_ElementsNumber * 3 + 1, _ElementsNumber * 3 + 1);
        var last = 0.0;
        
        
    
        Matrix<double> local = Matrix<double>.Build.Dense(4, 4);
        for (var i = 0; i < _ElementsNumber; i++)
        {
            
            double length = (double) _length / _ElementsNumber;
            
            if (_a != 0)
            {
                Matrix<double> aMat = DenseMatrix.OfArray(new[,]
                {
                    { 37 / (10 * length), -189 / (40 * length), 27 / (20 * length), -13 / (40 * length) },
                    { -189 / (40 * length), 54 / (5 * length), -297 / (40 * length), 27 / (20 * length) },
                    { 27 / (20 * length), -297 / (40 * length), 54 / (5 * length), -189 / (40 * length) },
                    { -13 / (40 * length), 27 / (20 * length), -189 / (40 * length), 37 / (10 * length) }
                });
                local = -_a * aMat;
            }
            
            
            
            if (_b != 0)
            {
                Matrix<double> bMat = DenseMatrix.OfArray(new[,]
                {
                    { -0.5, 57.0 / 80, -0.3, 7.0 / 80 },
                    { -57.0 / 80, 0, 81.0 / 80, -0.3 },
                    { 0.3, -81.0 / 80, 0, 57.0 / 80 },
                    { -7.0 / 80, 0.3, -57.0 / 80, 0.5 }
                });
                local += _b * bMat;
            }

            if (_c != 0)
            {
                Matrix<double> cMat = DenseMatrix.OfArray(new[,]
                {
                    { 8 * length / 105, 33 * length / 560, - 3 * length / 140 , 19 * length / 1680 },
                    { 33 * length / 560, 27 * length / 70, - 27 * length / 560, - 3 * length / 140 },
                    { - 3 * length / 140, -27 * length / 560, 27 * length / 70, 33 * length / 560 },
                    { 19 * length / 1680, - 3 * length / 140, 33 * length / 560, 8 * length / 105 }
                });
                local += _c * cMat;
            }
            
            
            
            if (i != 0)
            { 
                local[0, 0] += last;
                chm.SetSubMatrix(i * 3,i * 3, local);
            }
            else
                chm.SetSubMatrix(i, i, local);
            Debug.WriteLine(chm);
            last = local.Diagonal().Last();
        }
        

        return chm;
    }


    private Vector<double> CubForceMat()
    {


        Vector<double> cfm = Vector<double>.Build.Sparse(_ElementsNumber * 3 + 1);
        var last = 0.0;


        var k = 0;

        for (var i = 0; i < _ElementsNumber; i++)
        {
            //double l = Length(_mesh[k], _mesh[k + 1]) ;
            var length = (double) _length / _ElementsNumber;
            k += 1;




            Vector<double> local = Vector<double>.Build.Dense([
                _d * length / 8,
                _d * 3 * length / 8,
                _d * 3 * length / 8,
                _d * length / 8
            ]);
            
            
            
            
            if (i != 0)
            { 
                local[0] += last;
                cfm.SetSubVector(i * 3, 4, local);
            }
            else
                cfm.SetSubVector(i, 4, local);
            last = local.Last();
        }
        

        return cfm;
    }


    
        
    public double[] ErrorResult(Func<double, double> analyticalFunc, double[] numericalSol)
    {
        var mre = new double[_PointsNumber];
    
        for (int i = 0; i < _PointsNumber; i++)
            mre[i] = Math.Abs((analyticalFunc(_grid[i]) - numericalSol[i]) / analyticalFunc(_grid[i]));
    
        return mre;
    }



    public void ExportResult(Func<double, double> analyticalFunc, double[] linSol, double[] cubSol, double[] x, string name)
    {
        using var writer = new StreamWriter(Environment.GetFolderPath(Environment.SpecialFolder.UserProfile) 
                                            + "\\Downloads" + $"\\{name}", false, Encoding.UTF8);
        var n = _PointsNumber;
        writer.WriteLine($"x;" +
                         $"Аналитическое; " +
                         $"Линейное при {n - 1} КЭ; " +
                         $"Ошибка между Аналитическим и Линейным;" +
                         $"Кубическое при {n - 1} КЭ; " +
                         $"Ошибка между Аналитическим и Кубическим;" +
                         "Ошибка между формами");
        var errorl = ErrorResult(analyticalFunc, linSol);
        var errorc = ErrorResult(analyticalFunc, cubSol);
        
        for (var i = 0; i < n - 2; i++)
        {
            writer.WriteLine($"{x[i]};{analyticalFunc(x[i])};{linSol[i]};{errorl[i]};{cubSol[i]};{errorc[i]};" +
                             $"{Abs(errorl[i] - errorc[i])}");
        }
            
        
        
        
    }


    public double Length(double x1, double x2)
    {
        return Abs(x1 - x2);
    }


    public double[] GetMesh()
    {
        return _grid;
    }

    


    public int MinFinElemsWithSpecificRelativeError(double error, Func<double, double> func, double delta, int limit)
    {
        var prev = _ElementsNumber;
        for (int m = 20; m < limit; m++)
        {
            _ElementsNumber = m;
            _PointsNumber = m + 1;
            var sol = LinSolve();
            var rel = ErrorResult(func, sol.ToArray()).Max();
            if (Abs(rel - error) <= delta)
                return m;
        }

        Console.WriteLine("Out of bounds!");
        _ElementsNumber = prev;
        _PointsNumber = prev + 1;
        return 0;


    }


   

}



public class InitialConditions
{
    public int A { get; set; }
    public int B { get; set; }
    public int C { get; set; }
    public int D { get; set; }
    public int StartX { get; set; }
    public int EndX { get; set; }
    public EdgeCondition LeftCondition { get; set; }
    public EdgeCondition RightCondition { get; set; }
    
}
