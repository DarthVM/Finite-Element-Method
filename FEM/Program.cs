namespace FEM;

static class Program
{
    private static readonly InitialConditions Init = new InitialConditions
    {
        A = 32,
        B = 9,
        C = 0,
        D = 23,
        StartX = -2,
        EndX = 7,
        LeftCondition = new EdgeCondition
        {
            Side = ConditionSide.Left,
            Type = ConditionType.Third,
            Value = 32 // Такое же значение, что и в A для учета в матрице
        },
        RightCondition = new EdgeCondition
        {
            Side = ConditionSide.Right,
            Type = ConditionType.First,
            Value = -5
        }
    };


    [STAThread]
    static void Main()
    {
        Application.EnableVisualStyles();
        Application.SetCompatibleTextRenderingDefault(false);



        AnalyticalSolver analyticalSolver = new AnalyticalSolver(Init);
        var analyticalSol = analyticalSolver.Solve();
        
        FemSolver solver20 = new (Init, 20);

        var numLinSol20 = solver20.LinSolve().ToArray();
        var linMre20 = solver20.ErrorResult(analyticalSol, numLinSol20);
        var linMesh20 = solver20.GetMesh();
        
        
        var numCubSol20 = solver20.CubSolve().ToArray();
        var cubMre20 = solver20.ErrorResult(analyticalSol, numCubSol20);
        var cubMesh20 = solver20.GetMesh();
        
        Form1 form1 = new Form1(
            "Решение для 20 конечных элементов", 
            analyticalSol, 
            numLinSol20,
            linMesh20,
            numCubSol20,
            cubMesh20);
        
        
        

       
        
        
        FemSolver solver40 = new (Init, 40);

        var numLinSol40 = solver40.LinSolve().ToArray();
        var linMre40 = solver40.ErrorResult(analyticalSol, numLinSol40);
        var linMesh40 = solver40.GetMesh();

        
        var numCubSol40 = solver40.CubSolve().ToArray();
        var cubMre40 = solver40.ErrorResult(analyticalSol, numCubSol40);
        var cubMesh40 = solver40.GetMesh();
        
        
       
       
        Form1 form2 = new Form1(
            "Решение при 40 конечных элементов",
            analyticalSol,
            numLinSol40,
            linMesh40,
            numCubSol40,
            cubMesh40);


        Console.WriteLine($"Макс. ошибка линейной формы при 20 КЭ = {linMre20.Max()}");
        Console.WriteLine($"Макс. ошибка кубической формы при 20 КЭ = {cubMre20.Max()}");
        Console.WriteLine($"Макс. ошибка линейной формы при 20 КЭ = {linMre40.Max()}");
        Console.WriteLine($"Макс. ошибка кубической формы при 40 КЭ = {cubMre40.Max()}");
        
        
        
        solver20.ExportResult(analyticalSol, numLinSol20, numCubSol20, linMesh20, "Error20.csv");
        solver40.ExportResult(analyticalSol, numLinSol40, numCubSol40, linMesh40, "Error40.csv");
        
        
        Application.Run(new MultiFormContext(form1, form2));
    }
}




public class MultiFormContext : ApplicationContext
{
    private int openForms;
    public MultiFormContext(params Form[] forms)
    {
        openForms = forms.Length;

        foreach (var form in forms)
        {
            form.FormClosed += (s, args) =>
            {
                //When we have closed the last of the "starting" forms, 
                //end the program.
                if (Interlocked.Decrement(ref openForms) == 0)
                    ExitThread();
            };

            form.Show();
        }
    }
}