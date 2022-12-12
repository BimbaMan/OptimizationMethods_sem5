import java.util.Vector;

public abstract class Data
{
    public int N;

    public int getN(){
        return N;
    }
    public static boolean CheckConstraints(Vector<Double> vector) {
        return true;
    }
    public abstract double Function(Vector<Double> parameters);
}
