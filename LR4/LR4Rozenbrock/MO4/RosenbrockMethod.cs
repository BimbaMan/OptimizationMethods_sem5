using System;
using System.Collections.Generic;
using MathNet.Numerics.LinearAlgebra;

namespace RosenbrockMethod
{
   
    public sealed class RosenbrockMethod
    {
        private static Data input;

        private static int n;
        private static double beta;
        private static double alpha;
        private static double epsilon;
        private static int N;
        private static int l;

        private static Vector<double> delta0;

        public RosenbrockMethod(Data input)
        {
            n = input.N - 1;
            RosenbrockMethod.input = input;

            delta0 = Vector<double>.Build.Dense(input.N, 1.0); //начальная длина шага
            beta = -0.5; //коэф. сжатия
            alpha = 3; //коэф. растяжения
            epsilon = 0.1; //число для остановки алгоритма
            N = 1000;//количество итераций
        }
        
        public Vector<double> MinimizeFunction()
        {
            Matrix<double> d = Matrix<double>.Build.DenseDiagonal(n + 1, n + 1, 1.0);
            List<Vector<double>> x = new List<Vector<double>>{ Vector<double>.Build.Dense(n + 1, 30.0) };
            List<Vector<double>> y = new List<Vector<double>>{ x[0] };

            Vector<double> lambda = Vector<double>.Build.Dense(n + 1);
            Vector<double> delta = delta0.Clone();

            int iteration = 0;
            int currStep = 2;
            int i = 0;

            while (currStep != 0)
            {
                if (currStep == 2)
                {
                    currStep = Step2(ref delta, ref y, ref d, i);
                }
                if (currStep == 3)
                {
                    currStep = Step3(ref delta, ref y, ref x, ref i, iteration);
                }
                if (currStep == 4)
                {
                    currStep = Step4(ref delta, ref y, ref x, ref lambda, ref d, ref i, ref iteration);
                }
            }

            foreach (var item in x)
            {
                Console.WriteLine("\n" + item);
            }

            return x[^1];
        }

        private static int Step2(ref Vector<double> delta, ref List<Vector<double>> y, ref Matrix<double> d, int i) //шаг по iтому направлению 
        {
            Vector<double> temp = y[i].Add(d.Row(i).Multiply(delta[i])); //yi + deltai * di
            double result = input.Function(temp); //f(temp)
                             //f(yi)
            if (result < input.Function(y[i]) && input.CheckConstraints(temp))
            {
                l = 0;
                y.Add(temp); //yi+1 = temp
                delta[i] *= alpha; //delta = delta * alpha
            }
            else
            {
                l++;
                y.Add(y[i]); //yi+1 = yi
                delta[i] *= beta; //delta = delta * betta
            }

            return 3;
        }

        private static int Step3(ref Vector<double> delta, ref List<Vector<double>> y, ref List<Vector<double>> x, ref int i, int iteration)
        {
            if (i < n)
            {
                i++;
                return 2; //возврат к шагу 2
            }

            if (i == n) //проверить успешность поиска по текущим ортогональным направлениям
            {
                if (input.Function(y[n + 1]) < input.Function(y[0]) && input.CheckConstraints(y[n + 1])) 
                { // если хотя бы один шаг по направлению на шаге 2 был успешным, то 
                    Vector<double> temp = y[n + 1];//y1 = y(n+1)
                    y.Clear();
                    y.Add(temp); 
                    i = 0; //i =1 
                    return 2;
                }

                if (input.Function(y[n + 1]).Equals(input.Function(y[0]))) 
                { //если каждый из n последних шагов был неудачным оценить успешноть поиска по текущей итерации 
                    if (input.Function(y[n + 1]) < input.Function(x[iteration]) && input.CheckConstraints(y[n + 1]))
                    { //на kой итерации хотябы один шаг был удачный, то
                        return 4;// перейти к шагу 4
                    }

                    if (input.Function(y[n + 1]).Equals(input.Function(x[iteration])) && input.CheckConstraints(y[n + 1]))
                    {// не было ни одного удачного шага на кой итерации, процесс поиска приостановить 
                        if (l <= N) // если l последовательно неудачных серий шагов по всем направлениям на тек. итерации не превышает N, то
                        {
                            if (CheckEndCondition(ref delta)) //проверить условие окончания
                            {
                                return 0; //завершить поиск
                            }
                            else
                            {
                                Vector<double> temp = y[n + 1];
                                y.Clear();
                                y.Add(temp); //y1 = y(n+1)
                                i = 0;
                                return 2;
                            }
                        }
                    }
                }
            }

            return 4;
        }

        private static int Step4(ref Vector<double> delta, ref List<Vector<double>> y, ref List<Vector<double>> x, ref Vector<double> lambda, ref Matrix<double> d, ref int i, ref int iteration)
        {
            x.Add(y[n + 1]);
            iteration++;

            if (CheckEndCondition(ref x, ref lambda, iteration))
            {
                return 0; // завершить поиск
            }
            else
            {
                lambda = d.Solve(lambda); //вычислить джины шагов по каждому направлению поиска на какой итерации из соотношения x(k+1) - xk = СУММ(lamba*di)
                ProcedureGrammaShmita(ref lambda, ref d); //построить новый набор линейно независимых и взаино ортогональных направлений поиска d1...di c помощью процедуры Гамма-Шмидта
                delta = delta0.Clone();

                y.Clear();

                y.Add(x[iteration]);
                i = 0;
                l = 0;

                return 2; //переход к шагу 2 
            }
        }

        private static bool CheckEndCondition(ref List<Vector<double>> x, ref Vector<double> lambda, int iteration)
        {
            lambda = x[iteration].Subtract(x[iteration - 1]); 

            return lambda.L2Norm() <= epsilon;
        }

        private static bool CheckEndCondition(ref Vector<double> delta) 
        {// проверяются величины delta[i] использованные во время последней серии шагов 
            foreach (var item in delta)
            {
                if (Math.Abs(item) > epsilon)
                {
                    return false;
                }
            }

            return true;
        }
        
        private static void ProcedureGrammaShmita(ref Vector<double> lambda, ref Matrix<double> d)
        {
            List<Vector<double>> a = new List<Vector<double>>();
            List<Vector<double>> b = new List<Vector<double>>();

            for (int i = 0; i < n + 1; i++)
            {
                if (lambda[i].Equals(0))
                {
                    a.Add(d.Row(i));
                }
                else
                {
                    Vector<double> temp = Vector<double>.Build.Dense(n + 1, 0.0);

                    for (int j = i; j < n + 1; j++)
                    {
                        temp = temp.Add(d.Row(j).Multiply(lambda[j]));
                    }

                    a.Add(temp);
                }
            }

            for (int i = 0; i < n + 1; i++)
            {
                if (i != 0)
                {
                    Vector<double> temp = Vector<double>.Build.Sparse(n + 1);

                    for (int j = 0; j < i; j++)
                    {
                        temp = temp.Add(d.Row(j).Multiply(a[i].ConjugateDotProduct(d.Row(j))));
                    }

                    temp = a[i].Subtract(temp);
                    b.Add(temp);
                }
                else
                {
                    b.Add(a[i]);
                }

                d.SetRow(i,b[i].Divide(b[i].L2Norm()));
            }
        }
    }
}