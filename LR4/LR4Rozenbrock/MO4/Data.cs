using MathNet.Numerics.LinearAlgebra;

namespace RosenbrockMethod
{
   
    public abstract class Data
    {
        public abstract int N { get; }
        public virtual bool CheckConstraints(Vector<double> vector)
        {
            return true;
        }
        public abstract double Function(Vector<double> parameters);
    }
}
