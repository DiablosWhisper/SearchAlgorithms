using System.Collections.Generic;
using Console = Colorful.Console;
using System.Drawing;
using System.Linq;
using System;

namespace Graphs
{
    public class Utilities
    {
        protected double Infinity { get { return double.PositiveInfinity; } }
        protected double Min(double FirstValue, double SecondValue)
        {
            if (FirstValue < SecondValue)
            {
                return FirstValue;
            }

            return SecondValue;
        }
        protected double CalculateHeuristic(Node Start, Node Target)
        {
            return Math.Sqrt(Math.Pow(Start.Coordinate.X - Target.Coordinate.X, 2) + Math.Pow(Start.Coordinate.Y - Target.Coordinate.Y, 2));
        }
        protected void SetNodesWeightOfInfinity(List<Node> Nodes)
        {
            for (int Index = 0; Index < Nodes.Count; Index++)
            {
                Nodes[Index].Weight = Infinity;
            }
        }
        public void PrintResult(List<Node> Nodes, Node Start)
        {
            Console.Write($"{Start.Name} => ");

            Console.WriteLine(string.Join($"{Start.Name} => ", Nodes.Select(Node => Node.Name + $" = {Node.Weight}\n")));
        }
        public void PrintResult(List<Edge> Edges)
        {
            Console.WriteLine(string.Join(" => ", Edges.Select(Edge => Edge.Name)));
        }
        public void PrintResult(double[,] Matrix)
        {
            Console.Write("\t");

            for (int Index = 0; Index < Matrix.GetLength(0); Index++)
            {
                Console.Write($"[{Index.ToString()}]\t");
            }

            Console.WriteLine("\n");

            for (int FirstIndex = 0; FirstIndex < Matrix.GetLength(0); FirstIndex++)
            {
                Console.Write($"[{FirstIndex}]\t");

                for (int SecondIndex = 0; SecondIndex < Matrix.GetLength(0); SecondIndex++)
                {
                    Console.Write($"{Matrix[FirstIndex, SecondIndex]}\t");
                }

                Console.WriteLine("\n");
            }
        }
    }

    public class Coordinate
    {
        public Coordinate(int X, int Y)
        {
            this.X = X;

            this.Y = Y;
        }
        public int X { get; private set; }
        public int Y { get; private set; }
        public Coordinate()
        {
            X = 0;

            Y = 0;
        }
    }

    public class Node : ICloneable
    {
        public Node(int NumberOfInheritors, string Name, Coordinate Coordinate)
        {
            Inheritors = new Node[NumberOfInheritors + 1];

            this.Coordinate = Coordinate;

            this.Name = $"[{Name}]";

            Weight = 0.0;

            Index = 0;
        }
        public Node(int NumberOfInheritors, string Name)
        {
            Inheritors = new Node[NumberOfInheritors + 1];

            this.Name = $"[{Name}]";

            Weight = 0.0;

            Index = 0;
        }
        public Coordinate Coordinate { get; private set; }
        public double HeuristicPathWeight { get; set; }
        public double GainedPathWeight { get; set; }
        public override bool Equals(object Object)
        {
            Node Node = (Node)Object;

            return this.Name.Equals(Node.Name);
        }
        public double TotalPathWeight { get; set; }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        public Node[] Inheritors { get; set; }
        public double Weight { get; set; }
        public string Name { get; set; }
        public Node this[int Index]
        {
            get
            {
                if (Index >= 0 && Index < Inheritors.Length)
                {
                    return Inheritors[Index];
                }

                return null;
            }
            set
            {
                if (Index >= 0 && Index < Inheritors.Length)
                {
                    Inheritors[Index] = value;
                }
            }
        }
        public int Index { get; set; }
        public object Clone()
        {
            Node Node = new Node
            {
                Coordinate = this.Coordinate,

                Inheritors = this.Inheritors,

                Weight = this.Weight,

                Name = this.Name,

                Index = this.Index,
            };

            return Node;
        }
        public Node()
        {
            HeuristicPathWeight = 0.0;

            GainedPathWeight = 0.0;

            TotalPathWeight = 0.0;

            Name = string.Empty;

            Coordinate = null;

            Inheritors = null;

            Weight = 0.0;

            Index = 0;
        }
    }

    public class Edge : ICloneable
    {
        public Edge(double Weight, params Node[] Nodes)
        {
            Ends = new List<Node>();

            for (int Index = 0; Index < Nodes.Length; Index++)
            {
                Ends.Add(Nodes[Index]);
            }

            this.Weight = Weight;

            Ends[0][Nodes[1].Index] = Nodes[1];

            Name = $"{Nodes[0].Name}{Nodes[1].Name}";
        }
        public List<Node> Ends { get; private set; }
        public string Name { get; private set; }
        public double Weight { get; set; }
        public Node this[int Index]
        {
            get
            {
                if (Index >= 0 && Index <= 1)
                {
                    return Ends[Index];
                }

                return null;
            }
        }
        public int Index { get; set; }
        public object Clone()
        {
            List<Node> EndsOfEdge = new List<Node>()
            {
                (Node)this.Ends[0].Clone(),

                (Node)this.Ends[1].Clone()
            };

            Edge Edge = new Edge
            {
                Weight = this.Weight,

                Name = this.Name,

                Index = this.Index,

                Ends = EndsOfEdge
            };

            return Edge;
        }
        public Edge()
        {
            Ends = new List<Node>();

            Name = string.Empty;

            Weight = 0.0;

            Index = 0;
        }
    }

    public class Graph : ICloneable
    {
        public int NumberOfNodes { get { return SetOfNodes.Count; } }
        public int NumberOfEdges { get { return SetOfEdges.Count; } }
        public bool NegativeCycleChecker(List<Node> Nodes)
        {
            for (int Index = 0; Index < NumberOfEdges; Index++)
            {
                if (Nodes[SetOfEdges[Index][1].Index].Weight > Nodes[SetOfEdges[Index][0].Index].Weight + SetOfEdges[Index].Weight)
                {
                    return false;
                }
            }
            return true;
        }
        public bool NegativeCycleChecker(Node[] Nodes)
        {
            for (int Index = 0; Index < NumberOfEdges; Index++)
            {
                if (Nodes[SetOfEdges[Index][1].Index].Weight > Nodes[SetOfEdges[Index][0].Index].Weight + SetOfEdges[Index].Weight)
                {
                    return false;
                }
            }
            return true;
        }
        public List<Node> SetOfNodes { get; private set; }
        public List<Edge> SetOfEdges { get; set; }
        public Edge FindEdge(params Node[] Nodes)
        {
            for (int Index = 0; Index < NumberOfEdges; Index++)
            {
                if (SetOfEdges[Index][0].Equals(Nodes[0]) && SetOfEdges[Index][1].Equals(Nodes[1]))
                {
                    return SetOfEdges[Index];
                }
            }

            return new Edge { Weight = double.PositiveInfinity };
        }
        public void AddOneWayEdge(Edge Edge)
        {
            SetOfEdges.Add(Edge);

            SetOfEdges[SetOfEdges.Count - 1].Index = SetOfEdges.Count - 1;
        }
        public void AddTwoWayEdge(Edge Edge)
        {
            SetOfEdges.Add(Edge);

            SetOfEdges[SetOfEdges.Count - 1].Index = SetOfEdges.Count - 1;

            SetOfEdges.Add(new Edge(Edge.Weight, Edge[1], Edge[0]));

            SetOfEdges[SetOfEdges.Count - 1].Index = SetOfEdges.Count - 1;
        }
        public void AddNode(Node Node)
        {
            SetOfNodes.Add(Node);

            SetOfNodes[SetOfNodes.Count - 1].Index = SetOfNodes.Count - 1;
        }
        public object Clone()
        {
            List<Node> SetOfNodes = new List<Node>();

            for (int Index = 0; Index < NumberOfNodes; Index++)
            {
                SetOfNodes.Add((Node)this.SetOfNodes[Index].Clone());
            }

            List<Edge> SetOfEdges = new List<Edge>();

            for (int Index = 0; Index < NumberOfEdges; Index++)
            {
                SetOfEdges.Add((Edge)this.SetOfEdges[Index].Clone());
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
        public bool DepthFirstSearch(Node Start, Node Target, Node[] FoundPath)
        {
            VisitedNodes.Add(Start);

            if (VisitedNodes.Contains(Target))
            {
                return true;
            }

            for (int Index = 0; Index < InnerGraph.NumberOfNodes; Index++)
            {
                if (InnerGraph.SetOfNodes[Start.Index][Index] != null && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]) && InnerGraph.FindEdge(Start, InnerGraph.SetOfNodes[Index]).Weight > 0)
                {
                    FoundPath[Index] = Start;

                    if (DepthFirstSearch(InnerGraph.SetOfNodes[Index], Target, FoundPath))
                    {
                        return true;
                    }
                }
            }

            return VisitedNodes.Contains(Target);
        }
        public bool DepthFirstSearch(Node Start, Node Target)
        {
            VisitedNodes.Add(Start);

            if (VisitedNodes.Contains(Target))
            {
                return true;
            }

            for (int Index = 0; Index < InnerGraph.NumberOfNodes; Index++)
            {
                if (InnerGraph.SetOfNodes[Start.Index][Index] != null && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]))
                {
                    ResultOfSearching.Add(InnerGraph.FindEdge(Start, InnerGraph.SetOfNodes[Index]));

                    if (DepthFirstSearch(InnerGraph.SetOfNodes[Index], Target))
                    {
                        return true;
                    }
                }
            }

            return VisitedNodes.Contains(Target);
        }
        public DepthFirstSearchAlgorithm(Graph Graph)
        {
            InnerGraph = (Graph)Graph.Clone();

            ResultOfSearching = new List<Edge>();

            VisitedNodes = new List<Node>();
        }
        public bool DepthFirstSearch(Node Start)
        {
            VisitedNodes.Add(Start);

            for (int Index = 0; Index < InnerGraph.NumberOfNodes; Index++)
            {
                if (InnerGraph.SetOfNodes[Start.Index][Index] != null && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]))
                {
                    ResultOfSearching.Add(InnerGraph.FindEdge(Start, InnerGraph.SetOfNodes[Index]));

                    DepthFirstSearch(InnerGraph.SetOfNodes[Index]);
                }
            }

            return VisitedNodes.Count.Equals(InnerGraph.NumberOfNodes);
        }
        public List<Edge> ResultOfSearching { get; }
        private static Graph InnerGraph { get; set; }
        private List<Node> VisitedNodes { get; }
    }

    public class BreadthFirstSearchAlgorithm
    {
        public bool BreadthFirstSearch(Node Start, Node Target, Node[] FoundPath)
        {
            Queue.Enqueue(Start);

            VisitedNodes.Add(Start);

            FoundPath[Start.Index] = Start;

            if (VisitedNodes.Contains(Target))
            {
                return true;
            }

            for (; Queue.Count > 0;)
            {
                Start = Queue.Dequeue();

                for (int Index = 0; Index < InnerGraph.NumberOfNodes; Index++)
                {
                    if (InnerGraph.SetOfNodes[Start.Index][Index] != null && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]) && InnerGraph.FindEdge(Start, InnerGraph.SetOfNodes[Index]).Weight > 0)
                    {
                        FoundPath[Index] = Start;

                        Queue.Enqueue(InnerGraph.SetOfNodes[Index]);

                        VisitedNodes.Add(InnerGraph.SetOfNodes[Index]);
                    }
                }
            }

            return VisitedNodes.Contains(Target);
        }
        public bool BreadthFirstSearch(Node Start, Node Target)
        {
            Queue.Enqueue(Start);

            VisitedNodes.Add(Start);

            if (VisitedNodes.Contains(Target))
            {
                return true;
            }

            for (; Queue.Count > 0;)
            {
                Start = Queue.Dequeue();

                for (int Index = 0; Index < InnerGraph.NumberOfNodes; Index++)
                {
                    if (InnerGraph.SetOfNodes[Start.Index][Index] != null && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]))
                    {
                        ResultOfSearching.Add(InnerGraph.FindEdge(Start, InnerGraph.SetOfNodes[Index]));

                        Queue.Enqueue(InnerGraph.SetOfNodes[Index]);

                        VisitedNodes.Add(InnerGraph.SetOfNodes[Index]);

                        if (VisitedNodes.Contains(Target))
                        {
                            return true;
                        }
                    }
                }
            }

            return VisitedNodes.Contains(Target);
        }
        public BreadthFirstSearchAlgorithm(Graph Graph)
        {
            Queue = new Queue<Node>();

            InnerGraph = (Graph)Graph.Clone();

            ResultOfSearching = new List<Edge>();

            VisitedNodes = new List<Node>();
        }
        public bool BreadthFirstSearch(Node Start)
        {
            Queue.Enqueue(Start);

            VisitedNodes.Add(Start);

            for (; Queue.Count > 0;)
            {
                Start = Queue.Dequeue();

                for (int Index = 0; Index < InnerGraph.NumberOfNodes; Index++)
                {
                    if (InnerGraph.SetOfNodes[Start.Index][Index] != null && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]))
                    {
                        ResultOfSearching.Add(InnerGraph.FindEdge(Start, InnerGraph.SetOfNodes[Index]));

                        Queue.Enqueue(InnerGraph.SetOfNodes[Index]);

                        VisitedNodes.Add(InnerGraph.SetOfNodes[Index]);
                    }
                }
            }

            return VisitedNodes.Count.Equals(InnerGraph.NumberOfNodes);
        }
        public List<Edge> ResultOfSearching { get; }
        private static Graph InnerGraph { get; set; }
        private List<Node> VisitedNodes { get; }
        private Queue<Node> Queue { get; }
    }

    public class KruscalAlgorithm
    {
        public List<Edge> MinimumSpanningTree { get; }
        private static Graph InnerGraph { get; set; }
        public bool FindMinimumSpanningTree()
        {
            InnerGraph.SetOfEdges = InnerGraph.SetOfEdges.OrderBy(Edge => Edge.Weight).ToList<Edge>();

            for (int Index = 0; Index < InnerGraph.NumberOfEdges; Index++)
            {
                Node StartNode = FindSubTree(InnerGraph.SetOfEdges[Index][0]);

                Node EndNode = FindSubTree(InnerGraph.SetOfEdges[Index][1]);

                if (!StartNode.Equals(EndNode))
                {
                    MinimumSpanningTree.Add(InnerGraph.SetOfEdges[Index]);

                    InnerGraph.SetOfNodes[EndNode.Index] = StartNode;
                }
            }

            return true;
        }
        public KruscalAlgorithm(Graph Graph)
        {
            MinimumSpanningTree = new List<Edge>();

            InnerGraph = (Graph)Graph.Clone();
        }
        private Node FindSubTree(Node Root)
        {
            Node OldNode = new Node();

            Node RootOfTree = Root;

            for (; !RootOfTree.Equals(InnerGraph.SetOfNodes[RootOfTree.Index]);)
            {
                RootOfTree = InnerGraph.SetOfNodes[RootOfTree.Index];
            }

            for (; !Root.Equals(RootOfTree);)
            {
                OldNode = InnerGraph.SetOfNodes[Root.Index];

                InnerGraph.SetOfNodes[Root.Index] = RootOfTree;

                Root = OldNode;
            }

            return RootOfTree;
        }
    }

    public class PrimAlgorithm : Utilities
    {
        public List<Edge> MinimumSpanningTree { get; }
        private static Graph InnerGraph { get; set; }
        public bool FindMinimumSpanningTree()
        {
            VisitedNodes.Add(InnerGraph.SetOfNodes[0]);

            for (; VisitedNodes.Count < InnerGraph.NumberOfNodes;)
            {
                Node Start = new Node();

                Node End = new Node();

                Edge CurrentEdge = new Edge { Weight = Infinity };

                for (int Index = 0; Index < InnerGraph.NumberOfNodes; Index++)
                {
                    if (VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]))
                    {
                        foreach (Node Inheritor in InnerGraph.SetOfNodes[Index].Inheritors.Where(Node => Node != null && !VisitedNodes.Contains(Node)))
                        {
                            if (CurrentEdge.Weight > InnerGraph.FindEdge(InnerGraph.SetOfNodes[Index], Inheritor).Weight)
                            {
                                CurrentEdge = InnerGraph.FindEdge(InnerGraph.SetOfNodes[Index], Inheritor);

                                Start = CurrentEdge[0];

                                End = CurrentEdge[1];
                            }
                        }
                    }
                }

                VisitedNodes.Add(End);

                MinimumSpanningTree.Add(InnerGraph.FindEdge(Start, End));
            }

            return VisitedNodes.Count.Equals(InnerGraph.NumberOfNodes);
        }
        private List<Node> VisitedNodes { get; }
        public PrimAlgorithm(Graph Graph)
        {
            MinimumSpanningTree = new List<Edge>();

            InnerGraph = (Graph)Graph.Clone();

            VisitedNodes = new List<Node>();
        }
    }

    public class BellmanFordAlgorithm : Utilities
    {
        public List<Node> ShortestPathes { get; private set; }
        public bool FindTheShortestPathes(Node Start)
        {
            ShortestPathes[Start.Index].Weight = 0;

            for (int FirstIndex = 1; FirstIndex < InnerGraph.NumberOfNodes - 1; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < InnerGraph.NumberOfEdges; SecondIndex++)
                {
                    ShortestPathes[InnerGraph.SetOfEdges[SecondIndex][1].Index].Weight = Min(ShortestPathes[InnerGraph.SetOfEdges[SecondIndex][1].Index].Weight, ShortestPathes[InnerGraph.SetOfEdges[SecondIndex][0].Index].Weight + InnerGraph.SetOfEdges[SecondIndex].Weight);
                }
            }

            return InnerGraph.NegativeCycleChecker(ShortestPathes);
        }
        public BellmanFordAlgorithm(Graph Graph)
        {
            InnerGraph = (Graph)Graph.Clone();

            ShortestPathes = InnerGraph.SetOfNodes;

            SetNodesWeightOfInfinity(ShortestPathes);
        }
        private static Graph InnerGraph { get; set; }
    }

    public class DijkstraAlgorithm : Utilities
    {
        public List<Node> ShortestPathes { get; private set; }
        public bool FindTheShortestPathes(Node Start)
        {
            ShortestPathes[Start.Index].Weight = 0;

            for (int FirstIndex = 0; FirstIndex < InnerGraph.NumberOfNodes - 1; FirstIndex++)
            {
                Node CurrentNode = new Node { Weight = Infinity };

                for (int SecondIndex = 0; SecondIndex < InnerGraph.NumberOfNodes; SecondIndex++)
                {
                    if (!VisitedNodes.Contains(ShortestPathes[SecondIndex]) && ShortestPathes[SecondIndex].Weight <= CurrentNode.Weight)
                    {
                        CurrentNode = ShortestPathes[SecondIndex];

                        Index = ShortestPathes[SecondIndex].Index;
                    }
                }

                VisitedNodes.Add(ShortestPathes[Index]);

                for (int SecondIndex = 0; SecondIndex < InnerGraph.NumberOfNodes; SecondIndex++)
                {
                    if (!VisitedNodes.Contains(ShortestPathes[SecondIndex]) && InnerGraph.SetOfNodes[Index][SecondIndex] != null && !ShortestPathes[Index].Weight.Equals(Infinity) && ShortestPathes[SecondIndex].Weight > ShortestPathes[Index].Weight + InnerGraph.FindEdge(InnerGraph.SetOfNodes[Index], InnerGraph.SetOfNodes[SecondIndex]).Weight)
                    {
                        ShortestPathes[SecondIndex].Weight = ShortestPathes[Index].Weight + InnerGraph.FindEdge(InnerGraph.SetOfNodes[Index], InnerGraph.SetOfNodes[SecondIndex]).Weight;
                    }
                }
            }

            return InnerGraph.NegativeCycleChecker(ShortestPathes);
        }
        private static Graph InnerGraph { get; set; }
        public DijkstraAlgorithm(Graph Graph)
        {
            InnerGraph = (Graph)Graph.Clone();

            ShortestPathes = InnerGraph.SetOfNodes;

            VisitedNodes = new List<Node>();

            SetNodesWeightOfInfinity(ShortestPathes);
        }
        private List<Node> VisitedNodes { get; }
        private int Index { get; set; }
    }

    public class FloydWarshallAlgorithm : Utilities
    {
        public double[,] MatrixOfTheShortesPathes { get; private set; }
        public FloydWarshallAlgorithm(Graph Graph)
        {
            MatrixOfTheShortesPathes = new double[Graph.NumberOfNodes, Graph.NumberOfNodes];

            InnerGraph = (Graph)Graph.Clone();

            for (int FirstIndex = 0; FirstIndex < InnerGraph.NumberOfNodes; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < InnerGraph.NumberOfNodes; SecondIndex++)
                {
                    if (InnerGraph.SetOfNodes[FirstIndex][SecondIndex] != null)
                    {
                        MatrixOfTheShortesPathes[FirstIndex, SecondIndex] = InnerGraph.FindEdge(InnerGraph.SetOfNodes[FirstIndex], InnerGraph.SetOfNodes[SecondIndex]).Weight;
                    }
                    else
                    {
                        MatrixOfTheShortesPathes[FirstIndex, SecondIndex] = Infinity;
                    }
                }
            }
        }
        private static Graph InnerGraph { get; set; }
        public bool FindAllTheShortestPathes()
        {
            for (int Index = 0; Index < InnerGraph.NumberOfNodes; Index++)
            {
                MatrixOfTheShortesPathes[Index, Index] = 0;
            }

            for (int FirstIndex = 0; FirstIndex < InnerGraph.NumberOfNodes; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < InnerGraph.NumberOfNodes; SecondIndex++)
                {
                    for (int ThirdIndex = 0; ThirdIndex < InnerGraph.NumberOfNodes; ThirdIndex++)
                    {
                        MatrixOfTheShortesPathes[SecondIndex, ThirdIndex] = Min(MatrixOfTheShortesPathes[SecondIndex, ThirdIndex], MatrixOfTheShortesPathes[SecondIndex, FirstIndex] + MatrixOfTheShortesPathes[FirstIndex, ThirdIndex]);
                    }
                }
            }

            return true;
        }
    }

    public class JohnsonAlgorithm
    {
        private static BellmanFordAlgorithm BellmanFordPathSearch { get; set; }
        public double[,] MatrixOfTheShortesPathes { get; private set; }
        private static DijkstraAlgorithm DijkstraPathSearch { get; set; }
        private static Graph AuxiliaryGraph { get; set; }
        private static Graph InnerGraph { get; set; }
        public JohnsonAlgorithm(Graph Graph)
        {
            MatrixOfTheShortesPathes = new double[Graph.NumberOfNodes, Graph.NumberOfNodes];

            AuxiliaryGraph = (Graph)Graph.Clone();

            InnerGraph = (Graph)Graph.Clone();
        }
        public bool FindAllTheShortestPathes()
        {
            AuxiliaryGraph.AddNode(new Node(AuxiliaryGraph.NumberOfNodes, "Temporary"));

            for (int Index = 0; Index < AuxiliaryGraph.NumberOfNodes - 1; Index++)
            {
                AuxiliaryGraph.AddTwoWayEdge(new Edge(0, AuxiliaryGraph.SetOfNodes[AuxiliaryGraph.NumberOfNodes - 1], AuxiliaryGraph.SetOfNodes[Index]));
            }

            BellmanFordPathSearch = new BellmanFordAlgorithm(AuxiliaryGraph);

            if (BellmanFordPathSearch.FindTheShortestPathes(AuxiliaryGraph.SetOfNodes[AuxiliaryGraph.NumberOfNodes - 1]))
            {
                for (int FirstIndex = 0; FirstIndex < InnerGraph.NumberOfNodes; FirstIndex++)
                {
                    DijkstraPathSearch = new DijkstraAlgorithm(InnerGraph);

                    DijkstraPathSearch.FindTheShortestPathes(InnerGraph.SetOfNodes[FirstIndex]);

                    for (int SecondIndex = 0; SecondIndex < InnerGraph.NumberOfNodes; SecondIndex++)
                    {
                        MatrixOfTheShortesPathes[FirstIndex, SecondIndex] = DijkstraPathSearch.ShortestPathes[SecondIndex].Weight;
                    }
                }

                return true;
            }
            else
            {
                for (int Index = 0; Index < InnerGraph.NumberOfEdges; Index++)
                {
                    InnerGraph.SetOfEdges[Index].Weight = InnerGraph.SetOfEdges[Index].Weight + BellmanFordPathSearch.ShortestPathes[InnerGraph.SetOfEdges[Index][0].Index].Weight - BellmanFordPathSearch.ShortestPathes[InnerGraph.SetOfEdges[Index][1].Index].Weight;
                }

                for (int FirstIndex = 0; FirstIndex < InnerGraph.NumberOfNodes; FirstIndex++)
                {
                    DijkstraPathSearch = new DijkstraAlgorithm(InnerGraph);

                    DijkstraPathSearch.FindTheShortestPathes(InnerGraph.SetOfNodes[FirstIndex]);

                    for (int SecondIndex = 0; SecondIndex < InnerGraph.NumberOfNodes; SecondIndex++)
                    {
                        MatrixOfTheShortesPathes[FirstIndex, SecondIndex] = DijkstraPathSearch.ShortestPathes[SecondIndex].Weight;
                    }
                }
            }

            return true;
        }
    }

    public class FordFulkersonAlgorithm : Utilities
    {
        public bool FindMaximumFlowWithDFSRealization(Node Start, Node Target)
        {
            if (Start.Equals(Target))
            {
                CapacityOfFlow = 0.0;

                return true;
            }

            Node CurrentNode = new Node();

            for (; DepthFirstSearch.DepthFirstSearch(Start, Target, FoundPath);)
            {
                double CurrentStreamCapacity = Infinity;

                for (int Index = Target.Index; Index != Start.Index; Index = FoundPath[Index].Index)
                {
                    CurrentNode = FoundPath[Index];

                    CurrentStreamCapacity = Min(CurrentStreamCapacity, InnerGraph.FindEdge(CurrentNode, InnerGraph.SetOfNodes[Index]).Weight);
                }

                for (int Index = Target.Index; Index != Start.Index; Index = FoundPath[Index].Index)
                {
                    InnerGraph.FindEdge(CurrentNode, InnerGraph.SetOfNodes[Index]).Weight -= CurrentStreamCapacity;

                    InnerGraph.FindEdge(InnerGraph.SetOfNodes[Index], CurrentNode).Weight += CurrentStreamCapacity;
                }

                CapacityOfFlow += CurrentStreamCapacity;

                DepthFirstSearch = new DepthFirstSearchAlgorithm(InnerGraph);
            }

            return true;
        }
        public bool FindMaximumFlowWithBFSRealization(Node Start, Node Target)
        {
            if (Start.Equals(Target))
            {
                CapacityOfFlow = 0.0;

                return true;
            }

            Node CurrentNode = new Node();

            for (; BreadthFirstSearch.BreadthFirstSearch(Start, Target, FoundPath);)
            {
                double CurrentStreamCapacity = Infinity;

                for (int Index = Target.Index; Index != Start.Index; Index = FoundPath[Index].Index)
                {
                    CurrentNode = FoundPath[Index];

                    CurrentStreamCapacity = Min(CurrentStreamCapacity, InnerGraph.FindEdge(CurrentNode, InnerGraph.SetOfNodes[Index]).Weight);
                }

                for (int Index = Target.Index; Index != Start.Index; Index = FoundPath[Index].Index)
                {
                    InnerGraph.FindEdge(CurrentNode, InnerGraph.SetOfNodes[Index]).Weight -= CurrentStreamCapacity;

                    InnerGraph.FindEdge(InnerGraph.SetOfNodes[Index], CurrentNode).Weight += CurrentStreamCapacity;
                }

                CapacityOfFlow += CurrentStreamCapacity;

                BreadthFirstSearch = new BreadthFirstSearchAlgorithm(InnerGraph);
            }

            return true;
        }
        private static BreadthFirstSearchAlgorithm BreadthFirstSearch { get; set; }
        private static DepthFirstSearchAlgorithm DepthFirstSearch { get; set; }
        public double CapacityOfFlow { get; private set; }
        public FordFulkersonAlgorithm(Graph Graph)
        {
            InnerGraph = (Graph)Graph.Clone();

            BreadthFirstSearch = new BreadthFirstSearchAlgorithm(InnerGraph);

            DepthFirstSearch = new DepthFirstSearchAlgorithm(InnerGraph);

            FoundPath = new Node[InnerGraph.NumberOfNodes];

            CapacityOfFlow = 0.0;
        }
        private static Graph InnerGraph { get; set; }
        private Node[] FoundPath { get; set; }
    }

    public class AStarAlgorithm : Utilities
    {
        public bool FindTheShortesPath(Node Start, Node Target)
        {
            Node CurrentNode = Start;

            OpenList.Add(CurrentNode);

            for (; OpenList.Count > 0;)
            {
                CurrentNode = OpenList[0];

                if (CurrentNode.Equals(Target))
                {
                    VisitedNodes.Add(Target);

                    for (int Index = 0; Index < VisitedNodes.Count - 1; Index++)
                    {
                        ShortesPath.Add(InnerGraph.FindEdge(VisitedNodes[Index], VisitedNodes[Index + 1]));
                    }

                    return true;
                }

                OpenList.Remove(CurrentNode);

                ClosedList.Add(CurrentNode);

                foreach (Node Inheritor in CurrentNode.Inheritors.Where(Node => Node != null && Node.Index != Node.Inheritors.Length - 1))
                {
                    if (!ClosedList.Contains(Inheritor))
                    {
                        if (!OpenList.Contains(Inheritor))
                        {
                            Inheritor[Inheritor.Index] = CurrentNode;

                            Inheritor.HeuristicPathWeight = CalculateHeuristic(Inheritor, Target);

                            Inheritor.GainedPathWeight = InnerGraph.FindEdge(CurrentNode, Inheritor).Weight;

                            Inheritor.TotalPathWeight = Inheritor.GainedPathWeight + Inheritor.HeuristicPathWeight;

                            OpenList.Add(Inheritor);

                            OpenList = OpenList.OrderBy(Node => Node.TotalPathWeight).ToList<Node>();
                        }
                    }
                }

                VisitedNodes.Add(CurrentNode);
            }

            return true;
        }
        private static Graph InnerGraph { get; set; }
        private List<Node> ClosedList { get; set; }
        private List<Node> VisitedNodes { get; }
        private List<Node> OpenList { get; set; }
        public AStarAlgorithm(Graph Graph)
        {
            InnerGraph = (Graph)Graph.Clone();

            VisitedNodes = new List<Node>();

            ShortesPath = new List<Edge>();

            ClosedList = new List<Node>();

            OpenList = new List<Node>();
        }
        public List<Edge> ShortesPath { get; }
    }

    class Program
    {
        private static Graph Graph = new Graph();

        private static Graph NewGraph = new Graph();

        private static Utilities Utilities = new Utilities();
        static void Main(string[] args)
        {
            Graph.AddNode(new Node(7, "0", new Coordinate(0, 0)));
            Graph.AddNode(new Node(7, "1", new Coordinate(50, 100)));
            Graph.AddNode(new Node(7, "2", new Coordinate(100, 0)));
            Graph.AddNode(new Node(7, "3", new Coordinate(0, 200)));
            Graph.AddNode(new Node(7, "4", new Coordinate(250, 200)));
            Graph.AddNode(new Node(7, "5", new Coordinate(50, 250)));
            Graph.AddNode(new Node(7, "6", new Coordinate(150, 250)));

            Graph.AddTwoWayEdge(new Edge(7, Graph.SetOfNodes[0], Graph.SetOfNodes[1]));
            Graph.AddTwoWayEdge(new Edge(5, Graph.SetOfNodes[0], Graph.SetOfNodes[3]));
            Graph.AddTwoWayEdge(new Edge(8, Graph.SetOfNodes[1], Graph.SetOfNodes[2]));
            Graph.AddTwoWayEdge(new Edge(7, Graph.SetOfNodes[1], Graph.SetOfNodes[4]));
            Graph.AddTwoWayEdge(new Edge(9, Graph.SetOfNodes[1], Graph.SetOfNodes[3]));
            Graph.AddTwoWayEdge(new Edge(5, Graph.SetOfNodes[2], Graph.SetOfNodes[4]));
            Graph.AddTwoWayEdge(new Edge(15, Graph.SetOfNodes[3], Graph.SetOfNodes[4]));
            Graph.AddTwoWayEdge(new Edge(6, Graph.SetOfNodes[3], Graph.SetOfNodes[5]));
            Graph.AddTwoWayEdge(new Edge(8, Graph.SetOfNodes[5], Graph.SetOfNodes[4]));
            Graph.AddTwoWayEdge(new Edge(11, Graph.SetOfNodes[5], Graph.SetOfNodes[6]));
            Graph.AddTwoWayEdge(new Edge(9, Graph.SetOfNodes[6], Graph.SetOfNodes[4]));

            NewGraph.AddNode(new Node(7, "0", new Coordinate(0, 0)));
            NewGraph.AddNode(new Node(7, "1", new Coordinate(50, 100)));
            NewGraph.AddNode(new Node(7, "2", new Coordinate(100, 0)));
            NewGraph.AddNode(new Node(7, "3", new Coordinate(0, 200)));
            NewGraph.AddNode(new Node(7, "4", new Coordinate(500, 100)));
            NewGraph.AddNode(new Node(7, "5", new Coordinate(50, 250)));
            NewGraph.AddNode(new Node(7, "6", new Coordinate(150, 250)));

            NewGraph.AddOneWayEdge(new Edge(7, NewGraph.SetOfNodes[0], NewGraph.SetOfNodes[1]));
            NewGraph.AddOneWayEdge(new Edge(5, NewGraph.SetOfNodes[0], NewGraph.SetOfNodes[3]));
            NewGraph.AddOneWayEdge(new Edge(8, NewGraph.SetOfNodes[1], NewGraph.SetOfNodes[2]));
            NewGraph.AddOneWayEdge(new Edge(7, NewGraph.SetOfNodes[1], NewGraph.SetOfNodes[4]));
            NewGraph.AddOneWayEdge(new Edge(9, NewGraph.SetOfNodes[1], NewGraph.SetOfNodes[3]));
            NewGraph.AddOneWayEdge(new Edge(5, NewGraph.SetOfNodes[2], NewGraph.SetOfNodes[4]));
            NewGraph.AddOneWayEdge(new Edge(15, NewGraph.SetOfNodes[3], NewGraph.SetOfNodes[4]));
            NewGraph.AddOneWayEdge(new Edge(6, NewGraph.SetOfNodes[3], NewGraph.SetOfNodes[5]));
            NewGraph.AddOneWayEdge(new Edge(8, NewGraph.SetOfNodes[5], NewGraph.SetOfNodes[4]));
            NewGraph.AddOneWayEdge(new Edge(11, NewGraph.SetOfNodes[5], NewGraph.SetOfNodes[6]));
            NewGraph.AddOneWayEdge(new Edge(9, NewGraph.SetOfNodes[4], NewGraph.SetOfNodes[6]));

            DepthFirstSearchAlgorithm DFS = new DepthFirstSearchAlgorithm(Graph);
            BreadthFirstSearchAlgorithm BFS = new BreadthFirstSearchAlgorithm(Graph);
            KruscalAlgorithm KruskalTreeSearch = new KruscalAlgorithm(Graph);
            PrimAlgorithm PrimTreeSearch = new PrimAlgorithm(Graph);
            BellmanFordAlgorithm BellmanFordPathSearch = new BellmanFordAlgorithm(Graph);
            DijkstraAlgorithm DijkstraPathSearch = new DijkstraAlgorithm(Graph);
            FloydWarshallAlgorithm FloydWarshallPathSearch = new FloydWarshallAlgorithm(Graph);
            JohnsonAlgorithm JohnsonPathSearch = new JohnsonAlgorithm(Graph);
            AStarAlgorithm AStarPathSearch = new AStarAlgorithm(Graph);

            Console.WriteLine("[Current Graph]\n", Color.Green);

            Console.WriteLine("[0]        [2]");
            Console.WriteLine("|\\         /|");
            Console.WriteLine("| \\ 7    8/ |");
            Console.WriteLine("|  \\     /  |5");
            Console.WriteLine("|5  \\   /   |");
            Console.WriteLine("|    [1]    |");
            Console.WriteLine("| 9 /   \\7  |");
            Console.WriteLine("|  /     \\  |");
            Console.WriteLine("| /   15  \\ |");
            Console.WriteLine("[3]--------[4]");
            Console.WriteLine(" \\        /  \\ ");
            Console.WriteLine("  \\     8/    \\9");
            Console.WriteLine(" 6 \\    /      \\");
            Console.WriteLine("    \\  /        \\");
            Console.WriteLine("    [5]---------[6]");
            Console.WriteLine("           11\n");

            Console.WriteLine("[1 : Depth First Search]\n", Color.Green);

            DFS.DepthFirstSearch(Graph.SetOfNodes[0]);

            Utilities.PrintResult(DFS.ResultOfSearching);

            Console.WriteLine();

            Console.WriteLine("[2 : Breadth First Search]\n", Color.Green);

            BFS.BreadthFirstSearch(Graph.SetOfNodes[0]);

            Utilities.PrintResult(BFS.ResultOfSearching);

            Console.WriteLine();

            FordFulkersonAlgorithm FordFulkersonStreamSearch = new FordFulkersonAlgorithm(NewGraph);

            Console.WriteLine($"[3 : Kruskal Tree Search]\n", Color.Green);

            KruskalTreeSearch.FindMinimumSpanningTree();

            Utilities.PrintResult(KruskalTreeSearch.MinimumSpanningTree);

            Console.WriteLine();

            Console.WriteLine($"[4 : Prim Tree Search]\n", Color.Green);

            PrimTreeSearch.FindMinimumSpanningTree();

            Utilities.PrintResult(PrimTreeSearch.MinimumSpanningTree);

            Console.WriteLine();

            Console.WriteLine("[5 : Bellman-Ford Path Search]\n", Color.Green);

            BellmanFordPathSearch.FindTheShortestPathes(Graph.SetOfNodes[0]);

            Utilities.PrintResult(BellmanFordPathSearch.ShortestPathes, Graph.SetOfNodes[0]);

            Console.WriteLine("[6 : Dijkstra Path Search]\n", Color.Green);

            DijkstraPathSearch.FindTheShortestPathes(Graph.SetOfNodes[0]);

            Utilities.PrintResult(DijkstraPathSearch.ShortestPathes, Graph.SetOfNodes[0]);

            FloydWarshallPathSearch.FindAllTheShortestPathes();

            Console.WriteLine("[7 : Floyd-Warshall Path Search]\n", Color.Green);

            Utilities.PrintResult(FloydWarshallPathSearch.MatrixOfTheShortesPathes);

            Console.WriteLine();

            JohnsonPathSearch.FindAllTheShortestPathes();

            Console.WriteLine("[8 : Johnson Path Search]\n", Color.Green);

            Utilities.PrintResult(JohnsonPathSearch.MatrixOfTheShortesPathes);

            Console.WriteLine();

            Console.WriteLine("[9 : Ford-Fulkerson Stream Search]\n", Color.Green);

            FordFulkersonStreamSearch.FindMaximumFlowWithDFSRealization(NewGraph.SetOfNodes[0], NewGraph.SetOfNodes[6]);

            Console.WriteLine($"Maximum stream capacity is : {FordFulkersonStreamSearch.CapacityOfFlow}");

            Console.WriteLine();

            Console.WriteLine("[10 : A* Path Search]\n", Color.Green);

            AStarPathSearch.FindTheShortesPath(Graph.SetOfNodes[0], Graph.SetOfNodes[6]);

            Utilities.PrintResult(AStarPathSearch.ShortesPath);

            Console.ReadKey();
        }
    }
}