using MathNet.Numerics.LinearAlgebra;
using RosenbrockMethod;

namespace MO4
{
    class Program
    {
        public static void Main()
        {
            Console.WriteLine("Без ограничений");
            var unlimitData = new OptimizeWithoutLimit();
            var method1 = new RosenbrockMethod.RosenbrockMethod(unlimitData);
            var temp = method1.MinimizeFunction();
            
            Console.WriteLine("F = " + unlimitData.GetArea(temp));
            Console.WriteLine("L = " + unlimitData.Function(temp));

          /*  Console.WriteLine("\nС ограничениями");
            var limitData = new OptimizeWithLimit();
            var method2 = new RosenbrockMethod.RosenbrockMethod(limitData);
            temp = method2.MinimizeFunction();
            
            Console.WriteLine("F = " + limitData.GetArea(temp));
            Console.WriteLine("L = " + limitData.Function(temp));*/
        }
    }
}