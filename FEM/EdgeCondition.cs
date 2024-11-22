using MathNet.Numerics.LinearAlgebra;


namespace FEM;

public enum ConditionType
{
    First,
    Second,
    Third
}

public enum ConditionSide
{
    Left,
    Right
}




public class EdgeCondition
{
    public ConditionSide Side { get; set; }
    
    
    public ConditionType Type { get; set; }
    
    public double Value { get; set; }

    public void Init(ref Matrix<double> A, ref Vector<double> b)
    {
        switch (Side)
        {
            case ConditionSide.Left:
                InitLeftSide(ref A, ref b);
                break;
            case ConditionSide.Right:
                InitRightSide(ref A, ref b);
                break;
            default:
                throw new ArgumentOutOfRangeException();
        }
    }
    
    
    

    private void InitLeftSide(ref Matrix<double> A, ref Vector<double> b)
    {
        if (Type == ConditionType.First)
        {

            A.ClearRow(0);
            A[0, 0] = 1;
            b[0] = -Value;
        }
        else if (Type == ConditionType.Second)
        {
            b[0] += Value; // Value * _a
        }
        else
        {
             
            Vector<double> derivative = Vector<double>.Build.Dense(A.ColumnCount + 1);
            Vector<double> zeros = Vector<double>.Build.Dense(A.ColumnCount);
                
            A = A.InsertColumn(0, zeros);
                
            derivative[0] = 1;
            derivative[1] = -1;
            A = A.InsertRow(0, derivative);
            A[1, 0] = -Value; // -_a
            b = Vector<double>.Build.DenseOfEnumerable(new double[] { 0 }.Concat(b));  
        }
    }
    
    private void InitRightSide(ref Matrix<double> A, ref Vector<double> b)
    {
        if (Type == ConditionType.First)
        {

            A.ClearRow(A.RowCount - 1);
            A[A.RowCount - 1, A.ColumnCount - 1] = 1;
            b[^1] = - Value;
        }
        else if (Type == ConditionType.Second)
        {
            b[^1] += Value;
        }
        else
        {
             
            Vector<double> derivative = Vector<double>.Build.Dense(A.ColumnCount + 1);
            Vector<double> zeros = Vector<double>.Build.Dense(A.ColumnCount);
                
            A = A.InsertColumn(A.ColumnCount - 1, zeros);
                
            derivative[A.ColumnCount - 2] = 1;
            derivative[A.ColumnCount - 1] = -1;
            A = A.InsertRow(A.RowCount - 1, derivative);
            A[A.RowCount - 2, 0] = -Value;
            b = Vector<double>.Build.DenseOfEnumerable(new double[] { 0 }.Concat(b));
        }
    } 
}