﻿using System.Collections.Generic;
using System.Linq;
using System.IO;
using System;

namespace Graphs
{
    public class Node : ICloneable
    {
        public static bool operator ==(Node FirstNode, Node SecondNode)
        {
            if (FirstNode.NameOfNode.Equals(SecondNode.NameOfNode))
            {
                return true;
            }

            return false;
        }
        public static bool operator !=(Node FirstNode, Node SecondNode)
        {
            if (!FirstNode.NameOfNode.Equals(SecondNode.NameOfNode))
            {
                return true;
            }

            return false;
        }
        public override bool Equals(object Object)
        {
            Node Node = (Node)Object;

            return this.NameOfNode.Equals(Node.NameOfNode);
        }
        public double WeightOfNode { get; set; }
        public string NameOfNode { get; set; }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        public Node(int NumberOfNodes)
        {
            Inheritors = new bool[NumberOfNodes];

            NameOfNode = string.Empty;

            WeightOfNode = 0.0;

            IndexOfNode = 0;
        }
        public int IndexOfNode { get; set; }
        public bool[] Inheritors { get; set; }
        public object Clone()
        {
            Node Node = new Node
            {
                Inheritors = this.Inheritors,

                WeightOfNode = this.WeightOfNode,

                NameOfNode = this.NameOfNode,

                IndexOfNode = this.IndexOfNode
            };

            return Node;
        }
        public Node()
        {
            NameOfNode = string.Empty;

            Inheritors = null;

            WeightOfNode = 0.0;

            IndexOfNode = 0;
        }
    }

    public class Edge : ICloneable
    {
        public void CreateEdge(double WeightOfEdge, params Node[] Nodes)
        {
            for (int Index = 0; Index < Nodes.Length; Index++)
            {
                this.EndsOfEdge.Add(Nodes[Index]);
            }

            this.WeightOfEdge = WeightOfEdge;

            this.EndsOfEdge[0].Inheritors[Nodes[1].IndexOfNode] = true;

            this.EndsOfEdge[1].Inheritors[Nodes[0].IndexOfNode] = true;

            NameofEdge = Nodes[0].NameOfNode + Nodes[1].NameOfNode;
        }
        public List<Node> EndsOfEdge { get; private set; }
        public string NameofEdge { get; private set; }
        public double WeightOfEdge { get; set; }
        public int IndexOfEdge { get; set; }
        public Node this[int Index]
        {
            get
            {
                if (Index >= 0 && Index <= 1)
                {
                    return EndsOfEdge[Index];
                }

                return null;
            }
        }
        public object Clone()
        {
            List<Node> EndsOfEdge = new List<Node>()
            {
                (Node)this.EndsOfEdge[0].Clone(),

                (Node)this.EndsOfEdge[1].Clone()
            };

            Edge Edge = new Edge
            {
                WeightOfEdge = this.WeightOfEdge,

                NameofEdge = this.NameofEdge,

                IndexOfEdge = this.IndexOfEdge,

                EndsOfEdge = EndsOfEdge
            };

            return Edge;
        }
        public Edge()
        {
            EndsOfEdge = new List<Node>();

            NameofEdge = string.Empty;

            WeightOfEdge = 0.0;

            IndexOfEdge = 0;
        }
    }

    public class Graph : ICloneable
    {
        public Edge GetEdge(Node TheBeginningOfEdge, Node TheEndOfEdge)
        {
            foreach (Edge Edge in SetOfEdges.Where(CurrentEdge => CurrentEdge[0].Equals(TheBeginningOfEdge) && CurrentEdge[1].Equals(TheEndOfEdge) || CurrentEdge[0].Equals(TheEndOfEdge) && CurrentEdge[1].Equals(TheBeginningOfEdge)))
            {
                return Edge;
            }

            return new Edge
            {
                WeightOfEdge = double.PositiveInfinity
            };
        }
        public List<Node> SetOfNodes { get; set; }
        public List<Edge> SetOfEdges { get; set; }
        public object Clone()
        {
            List<Node> SetOfNodes = new List<Node>();

            foreach (var Node in this.SetOfNodes)
            {
                SetOfNodes.Add((Node)Node.Clone());
            }

            List<Edge> SetOfEdges = new List<Edge>();

            foreach (var Edge in this.SetOfEdges)
            {
                SetOfEdges.Add((Edge)Edge.Clone());
            }

            Graph Graph = new Graph
            {
                SetOfNodes = SetOfNodes,

                SetOfEdges = SetOfEdges
            };

            return Graph;
        }
        public Graph()
        {
            SetOfNodes = new List<Node>();

            SetOfEdges = new List<Edge>();
        }
    }

    public class DepthFirstSearchAlgorithm
    {
        public DepthFirstSearchAlgorithm(Graph Graph)
        {
            CopiedGraph = (Graph)Graph.Clone();

            VisitedNodes = new List<Node>();
        }
        public void DFS(Node TheBeginningOfSearch)
        {
            VisitedNodes.Add(TheBeginningOfSearch);

            for (int Index = 0; Index < CopiedGraph.SetOfNodes.Count; Index++)
            {
                if (CopiedGraph.SetOfNodes[TheBeginningOfSearch.IndexOfNode].Inheritors[Index].Equals(true) && !VisitedNodes.Contains(CopiedGraph.SetOfNodes[Index]))
                {
                    DFS(CopiedGraph.SetOfNodes[Index]);
                }
            }
        }
        private static Graph CopiedGraph { get; set; }
        public List<Node> VisitedNodes { get; }
    }

    public class BreadthFirstSearchAlgorithm
    {
        public BreadthFirstSearchAlgorithm(Graph Graph)
        {
            QueuedNodes = new Queue<Node>();

            CopiedGraph = (Graph)Graph.Clone();

            VisitedNodes = new List<Node>();
        }
        public void BFS(Node TheBeginningOfSearch)
        {
            QueuedNodes.Enqueue(TheBeginningOfSearch);

            VisitedNodes.Add(TheBeginningOfSearch);

            for (; QueuedNodes.Count > 0;)
            {
                TheBeginningOfSearch = QueuedNodes.Dequeue();

                for (int Index = 0; Index < CopiedGraph.SetOfNodes.Count; Index++)
                {
                    if (CopiedGraph.SetOfNodes[TheBeginningOfSearch.IndexOfNode].Inheritors[Index].Equals(true) && !VisitedNodes.Contains(CopiedGraph.SetOfNodes[Index]))
                    {
                        QueuedNodes.Enqueue(CopiedGraph.SetOfNodes[Index]);

                        VisitedNodes.Add(CopiedGraph.SetOfNodes[Index]);
                    }
                }
            }
        }
        private Queue<Node> QueuedNodes { get; }
        private static Graph CopiedGraph { get; set; }
        public List<Node> VisitedNodes { get; }
    }

    public class KruscalAlgorithm
    {
        public List<Edge> MinimumSpanningTree { get; }
        private static Graph CopiedGraph { get; set; }
        private Node FindSubTree(Node Node)
        {
            Node RootOfTree = Node;

            for (; !RootOfTree.Equals(CopiedGraph.SetOfNodes[RootOfTree.IndexOfNode]);)
            {
                RootOfTree = CopiedGraph.SetOfNodes[RootOfTree.IndexOfNode];
            }

            for (; !Node.Equals(RootOfTree);)
            {
                Node OldNode = CopiedGraph.SetOfNodes[Node.IndexOfNode];

                CopiedGraph.SetOfNodes[Node.IndexOfNode] = RootOfTree;

                Node = OldNode;
            }

            return RootOfTree;
        }
        public KruscalAlgorithm(Graph Graph)
        {
            MinimumSpanningTree = new List<Edge>();

            CopiedGraph = (Graph)Graph.Clone();
        }
        public void KruscalTreeSearch()
        {
            SortGraph();

            for (int Index = 0; Index < CopiedGraph.SetOfEdges.Count; Index++)
            {
                Node StartNode = FindSubTree(CopiedGraph.SetOfEdges[Index][0]);

                Node EndNode = FindSubTree(CopiedGraph.SetOfEdges[Index][1]);

                if (!StartNode.Equals(EndNode))
                {
                    MinimumSpanningTree.Add(CopiedGraph.SetOfEdges[Index]);

                    CopiedGraph.SetOfNodes[EndNode.IndexOfNode] = StartNode;
                }
            }
        }
        private void SortGraph()
        {
            Edge TemporaryEdge = new Edge();

            for (int FirstIndex = 0; FirstIndex < CopiedGraph.SetOfEdges.Count; FirstIndex++)
            {
                for (int SecondIndex = FirstIndex + 1; SecondIndex < CopiedGraph.SetOfEdges.Count; SecondIndex++)
                {
                    if (CopiedGraph.SetOfEdges[FirstIndex].WeightOfEdge > CopiedGraph.SetOfEdges[SecondIndex].WeightOfEdge)
                    {
                        TemporaryEdge = CopiedGraph.SetOfEdges[FirstIndex];

                        CopiedGraph.SetOfEdges[FirstIndex] = CopiedGraph.SetOfEdges[SecondIndex];

                        CopiedGraph.SetOfEdges[SecondIndex] = TemporaryEdge;
                    }
                }
            }
        }
    }

    public class PrimAlgorithm
    {
        public List<Node> VisitedNodes { get; private set; }
        public List<Edge> MinimumSpanningTree { get; }
        private static Graph CopiedGraph { get; set; }
        public PrimAlgorithm(Graph Graph)
        {
            MinimumSpanningTree = new List<Edge>();

            CopiedGraph = (Graph)Graph.Clone();

            VisitedNodes = new List<Node>();
        }
        public void PrimTreeSearch()
        {
            VisitedNodes.Add(CopiedGraph.SetOfNodes[0]);

            CopiedGraph.SetOfNodes.RemoveAt(VisitedNodes[0].IndexOfNode);

            for (; CopiedGraph.SetOfNodes.Count > 0;)
            {
                int IndexOfTheEasiestEdge = -1;

                for (int Index = 0; Index < CopiedGraph.SetOfEdges.Count; Index++)
                {
                    if (!VisitedNodes.IndexOf(CopiedGraph.SetOfEdges[Index][0]).Equals(-1) && !CopiedGraph.SetOfNodes.IndexOf(CopiedGraph.SetOfEdges[Index][1]).Equals(-1) || !VisitedNodes.IndexOf(CopiedGraph.SetOfEdges[Index][1]).Equals(-1) && !CopiedGraph.SetOfNodes.IndexOf(CopiedGraph.SetOfEdges[Index][0]).Equals(-1))
                    {
                        if (!IndexOfTheEasiestEdge.Equals(-1))
                        {
                            if (CopiedGraph.SetOfEdges[Index].WeightOfEdge < CopiedGraph.SetOfEdges[IndexOfTheEasiestEdge].WeightOfEdge)
                            {
                                IndexOfTheEasiestEdge = Index;
                            }
                        }
                        else
                        {
                            IndexOfTheEasiestEdge = Index;
                        }
                    }
                }

                if (!VisitedNodes.IndexOf(CopiedGraph.SetOfEdges[IndexOfTheEasiestEdge][0]).Equals(-1))
                {
                    VisitedNodes.Add(CopiedGraph.SetOfEdges[IndexOfTheEasiestEdge][1]);

                    CopiedGraph.SetOfNodes.Remove(CopiedGraph.SetOfEdges[IndexOfTheEasiestEdge][1]);
                }
                else
                {
                    VisitedNodes.Add(CopiedGraph.SetOfEdges[IndexOfTheEasiestEdge][0]);

                    CopiedGraph.SetOfNodes.Remove(CopiedGraph.SetOfEdges[IndexOfTheEasiestEdge][0]);
                }

                MinimumSpanningTree.Add(CopiedGraph.SetOfEdges[IndexOfTheEasiestEdge]);

                CopiedGraph.SetOfEdges.RemoveAt(IndexOfTheEasiestEdge);
            }
        }
    }

    public class BellmanFordAlgorithm
    {
        public void BellmanFordPathSearch(Node TheBeginningOfSearch)
        {
            foreach (Node Node in Pathes)
            {
                Pathes[Node.IndexOfNode].WeightOfNode = Infinity;
            }

            Pathes[TheBeginningOfSearch.IndexOfNode].WeightOfNode = 0;

            for (int FirstIndex = 1; FirstIndex < CopiedGraph.SetOfNodes.Count - 1; FirstIndex++)
            {
                foreach (var Edge in CopiedGraph.SetOfEdges.Where(CurrentEdge => Pathes[CopiedGraph.SetOfEdges[CurrentEdge.IndexOfEdge][1].IndexOfNode].WeightOfNode > Pathes[CopiedGraph.SetOfEdges[CurrentEdge.IndexOfEdge][0].IndexOfNode].WeightOfNode + CopiedGraph.SetOfEdges[CurrentEdge.IndexOfEdge].WeightOfEdge))
                {
                    Pathes[CopiedGraph.SetOfEdges[Edge.IndexOfEdge][1].IndexOfNode].WeightOfNode = Pathes[CopiedGraph.SetOfEdges[Edge.IndexOfEdge][0].IndexOfNode].WeightOfNode + CopiedGraph.SetOfEdges[Edge.IndexOfEdge].WeightOfEdge;
                }
            }
        }
        public BellmanFordAlgorithm(Graph Graph)
        {
            CopiedGraph = (Graph)Graph.Clone();

            Pathes = CopiedGraph.SetOfNodes;

            Infinity = double.PositiveInfinity;
        }
        public List<Node> Pathes { get; private set; }
        private static Graph CopiedGraph { get; set; }
        private double Infinity { get; set; }
    }

    public class DijkstraAlgorithm
    {
        public void DijkstraPathSearch(Node TheBeginningOfSearch)
        {
            foreach (Node Node in CopiedGraph.SetOfNodes)
            {
                Pathes[Node.IndexOfNode].WeightOfNode = Infinity;
            }

            Pathes[TheBeginningOfSearch.IndexOfNode].WeightOfNode = 0;

            for (int FirstIndex = 0; FirstIndex < CopiedGraph.SetOfNodes.Count - 1; FirstIndex++)
            {
                Node TemporaryNode = new Node
                {
                    WeightOfNode = Infinity
                };

                TemporaryNode.WeightOfNode = Infinity;

                for (int SecondIndex = 0; SecondIndex < CopiedGraph.SetOfNodes.Count; SecondIndex++)
                {
                    if (!VisitedNodes.Contains(Pathes[SecondIndex]) && Pathes[SecondIndex].WeightOfNode <= TemporaryNode.WeightOfNode)
                    {
                        TemporaryNode = Pathes[SecondIndex];

                        Index = Pathes[SecondIndex].IndexOfNode;
                    }
                }

                VisitedNodes.Add(Pathes[Index]);

                for (int SecondIndex = 0; SecondIndex < CopiedGraph.SetOfNodes.Count; SecondIndex++)
                {
                    if (!VisitedNodes.Contains(Pathes[SecondIndex]) && CopiedGraph.SetOfNodes[Index].Inheritors[SecondIndex].Equals(true) && !Pathes[Index].WeightOfNode.Equals(Infinity) && Pathes[SecondIndex].WeightOfNode > Pathes[Index].WeightOfNode + CopiedGraph.GetEdge(CopiedGraph.SetOfNodes[Index], CopiedGraph.SetOfNodes[SecondIndex]).WeightOfEdge)
                    {
                        Pathes[SecondIndex].WeightOfNode = Pathes[Index].WeightOfNode + CopiedGraph.GetEdge(CopiedGraph.SetOfNodes[Index], CopiedGraph.SetOfNodes[SecondIndex]).WeightOfEdge;
                    }
                }
            }
        }
        private static Graph CopiedGraph { get; set; }
        public List<Node> Pathes { get; private set; }
        private List<Node> VisitedNodes { get; }
        public DijkstraAlgorithm(Graph Graph)
        {
            CopiedGraph = (Graph)Graph.Clone();

            Pathes = CopiedGraph.SetOfNodes;

            VisitedNodes = new List<Node>();

            Infinity = double.PositiveInfinity;
        }
        private double Infinity { get; set; }
        private int Index { get; set; }
    }

    public class FloydWarshallAlgorithm
    {
        public double[,] MatrixOfShortesPathes { get; private set; }
        public FloydWarshallAlgorithm(Graph Graph)
        {
            MatrixOfShortesPathes = new double[Graph.SetOfNodes.Count, Graph.SetOfNodes.Count];

            CopiedGraph = (Graph)Graph.Clone();

            Infinity = double.PositiveInfinity;

            for (int FirstIndex = 0; FirstIndex < CopiedGraph.SetOfNodes.Count; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < CopiedGraph.SetOfNodes.Count; SecondIndex++)
                {
                    if (CopiedGraph.SetOfNodes[FirstIndex].Inheritors[SecondIndex].Equals(true))
                    {
                        MatrixOfShortesPathes[FirstIndex, SecondIndex] = CopiedGraph.GetEdge(CopiedGraph.SetOfNodes[FirstIndex], CopiedGraph.SetOfNodes[SecondIndex]).WeightOfEdge;
                    }
                    else
                    {
                        MatrixOfShortesPathes[FirstIndex, SecondIndex] = Infinity;
                    }
                }
            }
        }
        private static Graph CopiedGraph { get; set; }
        public void FloydWarshallPathSearch()
        {
            for (int Index = 0; Index < CopiedGraph.SetOfNodes.Count; Index++)
            {
                MatrixOfShortesPathes[Index, Index] = 0;
            }

            for (int FirstIndex = 0; FirstIndex < CopiedGraph.SetOfNodes.Count; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < CopiedGraph.SetOfNodes.Count; SecondIndex++)
                {
                    for (int ThirdIndex = 0; ThirdIndex < CopiedGraph.SetOfNodes.Count; ThirdIndex++)
                    {
                        if (MatrixOfShortesPathes[SecondIndex, ThirdIndex] > MatrixOfShortesPathes[SecondIndex, FirstIndex] + MatrixOfShortesPathes[FirstIndex, ThirdIndex])
                        {
                            MatrixOfShortesPathes[SecondIndex, ThirdIndex] = MatrixOfShortesPathes[SecondIndex, FirstIndex] + MatrixOfShortesPathes[FirstIndex, ThirdIndex];
                        }
                    }
                }
            }
        }
        private double Infinity { get; set; }
    }

    class Program
    {
        private static Graph Graph = new Graph();
        static void Main(string[] args)
        {
            for (int Index = 0; Index < 7; Index++)
            {
                Graph.SetOfNodes.Add(new Node(7));

                Graph.SetOfNodes[Index].IndexOfNode = Index;
            }

            for (int Index = 0; Index < 11; Index++)
            {
                Graph.SetOfEdges.Add(new Edge());

                Graph.SetOfEdges[Index].IndexOfEdge = Index;
            }

            Graph.SetOfNodes[0].NameOfNode = "A";
            Graph.SetOfNodes[1].NameOfNode = "B";
            Graph.SetOfNodes[2].NameOfNode = "C";
            Graph.SetOfNodes[3].NameOfNode = "D";
            Graph.SetOfNodes[4].NameOfNode = "E";
            Graph.SetOfNodes[5].NameOfNode = "F";
            Graph.SetOfNodes[6].NameOfNode = "G";

            Graph.SetOfEdges[0].CreateEdge(0, Graph.SetOfNodes[0], Graph.SetOfNodes[1]);
            Graph.SetOfEdges[1].CreateEdge(0, Graph.SetOfNodes[0], Graph.SetOfNodes[3]);
            Graph.SetOfEdges[2].CreateEdge(0, Graph.SetOfNodes[1], Graph.SetOfNodes[2]);
            Graph.SetOfEdges[3].CreateEdge(0, Graph.SetOfNodes[1], Graph.SetOfNodes[4]);
            Graph.SetOfEdges[4].CreateEdge(0, Graph.SetOfNodes[1], Graph.SetOfNodes[3]);
            Graph.SetOfEdges[5].CreateEdge(0, Graph.SetOfNodes[2], Graph.SetOfNodes[4]);
            Graph.SetOfEdges[6].CreateEdge(0, Graph.SetOfNodes[3], Graph.SetOfNodes[4]);
            Graph.SetOfEdges[7].CreateEdge(0, Graph.SetOfNodes[3], Graph.SetOfNodes[5]);
            Graph.SetOfEdges[8].CreateEdge(0, Graph.SetOfNodes[5], Graph.SetOfNodes[4]);
            Graph.SetOfEdges[9].CreateEdge(0, Graph.SetOfNodes[5], Graph.SetOfNodes[6]);
            Graph.SetOfEdges[10].CreateEdge(0, Graph.SetOfNodes[6], Graph.SetOfNodes[4]);

            Console.WriteLine("Current Graph\n");

            StreamReader Reader = new StreamReader("C:\\Users\\Gleb\\Desktop\\Graph.txt");

            Console.WriteLine(Reader.ReadToEnd());

            Reader.Close();

            DepthFirstSearchAlgorithm DFS = new DepthFirstSearchAlgorithm(Graph);

            DFS.DFS(Graph.SetOfNodes[0]);

            Console.WriteLine("1 : Depth First Search\n");

            Console.WriteLine(string.Join(" => ", DFS.VisitedNodes.Select(Vertex => Vertex.NameOfNode)));

            Console.WriteLine();

            BreadthFirstSearchAlgorithm BFS = new BreadthFirstSearchAlgorithm(Graph);

            BFS.BFS(Graph.SetOfNodes[0]);

            Console.WriteLine("2 : Breadth First Search\n");

            Console.WriteLine(string.Join(" => ", BFS.VisitedNodes.Select(Vertex => Vertex.NameOfNode)));

            Console.WriteLine();

            Graph.SetOfEdges[0].WeightOfEdge = 7;
            Graph.SetOfEdges[1].WeightOfEdge = 5;
            Graph.SetOfEdges[2].WeightOfEdge = 8;
            Graph.SetOfEdges[3].WeightOfEdge = 7;
            Graph.SetOfEdges[4].WeightOfEdge = 9;
            Graph.SetOfEdges[5].WeightOfEdge = 5;
            Graph.SetOfEdges[6].WeightOfEdge = 15;
            Graph.SetOfEdges[7].WeightOfEdge = 6;
            Graph.SetOfEdges[8].WeightOfEdge = 8;
            Graph.SetOfEdges[9].WeightOfEdge = 11;
            Graph.SetOfEdges[10].WeightOfEdge = 9;

            KruscalAlgorithm KruskalTreeSearch = new KruscalAlgorithm(Graph);

            KruskalTreeSearch.KruscalTreeSearch();

            Console.WriteLine("3 : Kruskal Tree Search\n");

            Console.WriteLine(string.Join(" => ", KruskalTreeSearch.MinimumSpanningTree.Select(Edges => Edges.NameofEdge)));

            Console.WriteLine();

            PrimAlgorithm PrimTreeSearch = new PrimAlgorithm(Graph);

            PrimTreeSearch.PrimTreeSearch();

            Console.WriteLine("4 : Prim Tree Search\n");

            Console.WriteLine(string.Join(" => ", PrimTreeSearch.MinimumSpanningTree.Select(Edges => Edges.NameofEdge)));

            Console.WriteLine();

            Console.WriteLine("5 : Bellman-Ford Path Search\n");

            BellmanFordAlgorithm BellmanFordPathSearch = new BellmanFordAlgorithm(Graph);

            BellmanFordPathSearch.BellmanFordPathSearch(Graph.SetOfNodes[0]);

            Console.Write($"{Graph.SetOfNodes[0].NameOfNode} => ");

            Console.WriteLine(string.Join($"{Graph.SetOfNodes[0].NameOfNode} => ", BellmanFordPathSearch.Pathes.Select(Nodes => Nodes.NameOfNode + $" = {Nodes.WeightOfNode}\n")));

            DijkstraAlgorithm DijkstraPathSearch = new DijkstraAlgorithm(Graph);

            DijkstraPathSearch.DijkstraPathSearch(Graph.SetOfNodes[0]);

            Console.WriteLine("6 : Dijkstra Path Search\n");

            Console.Write($"{Graph.SetOfNodes[0].NameOfNode} => ");

            Console.WriteLine(string.Join($"{Graph.SetOfNodes[0].NameOfNode} => ", DijkstraPathSearch.Pathes.Select(Nodes => Nodes.NameOfNode + $" = {Nodes.WeightOfNode}\n")));

            FloydWarshallAlgorithm FloydWarshallPathSearch = new FloydWarshallAlgorithm(Graph);

            FloydWarshallPathSearch.FloydWarshallPathSearch();

            Console.WriteLine("7 :Floyd-Warshall Path Search\n");

            Console.Write("\t");

            for (int Index = 0; Index < Graph.SetOfNodes.Count; Index++)
            {
                Console.Write($"{Graph.SetOfNodes[Index].NameOfNode}\t");
            }

            Console.WriteLine();

            for (int FirstIndex = 0; FirstIndex < Graph.SetOfNodes.Count; FirstIndex++)
            {
                Console.Write($"{Graph.SetOfNodes[FirstIndex].NameOfNode}\t");

                for (int SecondIndex = 0; SecondIndex < Graph.SetOfNodes.Count; SecondIndex++)
                {
                    Console.Write($"{FloydWarshallPathSearch.MatrixOfShortesPathes[FirstIndex, SecondIndex]} \t");
                }

                Console.WriteLine();
            }

            Console.ReadKey();
        }
    }
}