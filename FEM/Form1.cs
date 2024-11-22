using OxyPlot;
using OxyPlot.Series;
using OxyPlot.Axes;

namespace FEM;

public partial class Form1 : Form
{
    public Form1(
        string title, 
        Func<double, double> func, 
        double[] lin,
        double[] linX,
        double[] cub,
        double[] cubX
        )
    {
        InitializeComponent();
        var myModel = new PlotModel { Title = title };
        
        AddScatter(myModel, lin.ToList(), linX.ToList(), "Линейная форма", OxyColors.Automatic);
        AddScatter(myModel, cub.ToList(), cubX.ToList(), "Кубическая форма", OxyColors.Automatic);
        myModel.Series.Add(new FunctionSeries(func, linX[0], linX[^1], 1e-2, "Аналитическое решение"));
        
        
        
        plot1.Model = myModel;
    }

    public void AddScatter(PlotModel model, List<double> y, List<double> x, string name, OxyColor color)
    {
        var lineSeries = new LineSeries
        {
            Title = name,
            MarkerType = MarkerType.Diamond,
            MarkerSize = 4,
            MarkerStroke = color,
            MarkerFill = OxyColors.Automatic
        };

        model.Axes.Clear();
        
        model.Axes.Add(new LinearAxis
        {
            Position = AxisPosition.Left,
            Title = "U(x)"
        });
        
        model.Axes.Add(new LinearAxis
        {
            Position = AxisPosition.Bottom,
            Title = "X"
        });
        
        
        
        
        
        var length = x.Count;

        for (int i = 0; i < length; i++)
        {
            lineSeries.Points.Add(new DataPoint(x[i], y[i]));
        }
        
        
        model.Series.Add(lineSeries);
    }
}