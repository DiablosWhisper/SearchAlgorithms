using System.Collections.Generic;
using System.Linq;
using System;

namespace Graphs
{
    public class Node : ICloneable
    {
        public Node(int NumberOfInheritors, string Name)
        {
            Connections = new bool[NumberOfInheritors + 1];

            Weight = double.PositiveInfinity;

            this.Name = $"[{Name}]";

            Index = 0;
        }
        public double HeuristicPathWeight { get; set; }
        public double GainedPathWeight { get; set; }
        public override bool Equals(object Object)
        {
            Node Node = (Node)Object;

            return this.Name.Equals(Node.Name);
        }
        public double FinalPathWeight { get; set; }
        public bool[] Connections { get; set; }
        public override int GetHashCode()
        {
            return base.GetHashCode();
        }
        public double Weight { get; set; }
        public string Name { get; set; }
        public int Index { get; set; }
        public object Clone()
        {
            Node Node = new Node
            {
                Connections = this.Connections,

                Weight = this.Weight,

                Name = this.Name,

                Index = this.Index
            };

            return Node;
        }
        public Node()
        {
            Weight = double.PositiveInfinity;

            HeuristicPathWeight = 0.0;

            GainedPathWeight = 0.0;

            FinalPathWeight = 0.0;

            Name = string.Empty;

            Connections = null;

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

            Ends[0].Connections[Nodes[1].Index] = true;

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
        public List<Edge> SetOfEdges { get; private set; }
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
        private double Infinity { get; set; }
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

            Infinity = double.PositiveInfinity;
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
                if (InnerGraph.SetOfNodes[Start.Index].Connections[Index].Equals(true) && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]) && InnerGraph.FindEdge(Start, InnerGraph.SetOfNodes[Index]).Weight > 0)
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
                if (InnerGraph.SetOfNodes[Start.Index].Connections[Index].Equals(true) && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]))
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
                if (InnerGraph.SetOfNodes[Start.Index].Connections[Index].Equals(true) && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]))
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
                    if (InnerGraph.SetOfNodes[Start.Index].Connections[Index].Equals(true) && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]) && InnerGraph.FindEdge(Start, InnerGraph.SetOfNodes[Index]).Weight > 0)
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
                    if (InnerGraph.SetOfNodes[Start.Index].Connections[Index].Equals(true) && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]))
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
                    if (InnerGraph.SetOfNodes[Start.Index].Connections[Index].Equals(true) && !VisitedNodes.Contains(InnerGraph.SetOfNodes[Index]))
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
            SortGraph();

            Node StartNode = new Node();

            Node EndNode = new Node();

            for (int Index = 0; Index < InnerGraph.NumberOfEdges; Index++)
            {
                StartNode = FindSubTree(InnerGraph.SetOfEdges[Index][0]);

                EndNode = FindSubTree(InnerGraph.SetOfEdges[Index][1]);

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
        private void SortGraph()
        {
            Edge TemporaryEdge = new Edge();

            for (int FirstIndex = 0; FirstIndex < InnerGraph.NumberOfEdges; FirstIndex++)
            {
                for (int SecondIndex = FirstIndex + 1; SecondIndex < InnerGraph.NumberOfEdges; SecondIndex++)
                {
                    if (InnerGraph.SetOfEdges[FirstIndex].Weight > InnerGraph.SetOfEdges[SecondIndex].Weight)
                    {
                        TemporaryEdge = InnerGraph.SetOfEdges[FirstIndex];

                        InnerGraph.SetOfEdges[FirstIndex] = InnerGraph.SetOfEdges[SecondIndex];

                        InnerGraph.SetOfEdges[SecondIndex] = TemporaryEdge;
                    }
                }
            }
        }
    }

    public class PrimAlgorithm
    {
        public List<Edge> MinimumSpanningTree { get; }
        private static Graph InnerGraph { get; set; }
        public bool FindMinimumSpanningTree()
        {
            VisitedNodes.Add(InnerGraph.SetOfNodes[0]);

            InnerGraph.SetOfNodes.RemoveAt(VisitedNodes[0].Index);

            for (; InnerGraph.NumberOfNodes > 0;)
            {
                int IndexOfTheEasiestEdge = -1;

                for (int Index = 0; Index < InnerGraph.NumberOfEdges; Index++)
                {
                    if (!VisitedNodes.IndexOf(InnerGraph.SetOfEdges[Index][0]).Equals(-1) && !InnerGraph.SetOfNodes.IndexOf(InnerGraph.SetOfEdges[Index][1]).Equals(-1) || !VisitedNodes.IndexOf(InnerGraph.SetOfEdges[Index][1]).Equals(-1) && !InnerGraph.SetOfNodes.IndexOf(InnerGraph.SetOfEdges[Index][0]).Equals(-1))
                    {
                        if (!IndexOfTheEasiestEdge.Equals(-1))
                        {
                            if (InnerGraph.SetOfEdges[Index].Weight < InnerGraph.SetOfEdges[IndexOfTheEasiestEdge].Weight)
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

                if (!VisitedNodes.IndexOf(InnerGraph.SetOfEdges[IndexOfTheEasiestEdge][0]).Equals(-1))
                {
                    VisitedNodes.Add(InnerGraph.SetOfEdges[IndexOfTheEasiestEdge][1]);

                    InnerGraph.SetOfNodes.Remove(InnerGraph.SetOfEdges[IndexOfTheEasiestEdge][1]);
                }
                else
                {
                    VisitedNodes.Add(InnerGraph.SetOfEdges[IndexOfTheEasiestEdge][0]);

                    InnerGraph.SetOfNodes.Remove(InnerGraph.SetOfEdges[IndexOfTheEasiestEdge][0]);
                }

                MinimumSpanningTree.Add(InnerGraph.SetOfEdges[IndexOfTheEasiestEdge]);

                InnerGraph.SetOfEdges.RemoveAt(IndexOfTheEasiestEdge);
            }

            return true;
        }
        private List<Node> VisitedNodes { get; }
        public PrimAlgorithm(Graph Graph)
        {
            MinimumSpanningTree = new List<Edge>();

            InnerGraph = (Graph)Graph.Clone();

            VisitedNodes = new List<Node>();
        }
    }

    public class BellmanFordAlgorithm
    {
        public List<Node> ShortestPathes { get; private set; }
        public bool FindTheShortestPathes(Node Start)
        {
            ShortestPathes[Start.Index].Weight = 0;

            for (int FirstIndex = 1; FirstIndex < InnerGraph.NumberOfNodes - 1; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < InnerGraph.NumberOfEdges; SecondIndex++)
                {
                    ShortestPathes[InnerGraph.SetOfEdges[SecondIndex][1].Index].Weight = Math.Min(ShortestPathes[InnerGraph.SetOfEdges[SecondIndex][1].Index].Weight, ShortestPathes[InnerGraph.SetOfEdges[SecondIndex][0].Index].Weight + InnerGraph.SetOfEdges[SecondIndex].Weight);
                }
            }

            return InnerGraph.NegativeCycleChecker(ShortestPathes);
        }
        public BellmanFordAlgorithm(Graph Graph)
        {
            InnerGraph = (Graph)Graph.Clone();

            ShortestPathes = InnerGraph.SetOfNodes;

            Infinity = double.PositiveInfinity;
        }
        private static Graph InnerGraph { get; set; }
        private double Infinity { get; set; }
    }

    public class DijkstraAlgorithm
    {
        public List<Node> ShortestPathes { get; private set; }
        public bool FindTheShortestPathes(Node Start)
        {
            ShortestPathes[Start.Index].Weight = 0;

            for (int FirstIndex = 0; FirstIndex < InnerGraph.NumberOfNodes - 1; FirstIndex++)
            {
                Node CurrentNode = new Node();

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
                    if (!VisitedNodes.Contains(ShortestPathes[SecondIndex]) && InnerGraph.SetOfNodes[Index].Connections[SecondIndex].Equals(true) && !ShortestPathes[Index].Weight.Equals(Infinity) && ShortestPathes[SecondIndex].Weight > ShortestPathes[Index].Weight + InnerGraph.FindEdge(InnerGraph.SetOfNodes[Index], InnerGraph.SetOfNodes[SecondIndex]).Weight)
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

            Infinity = double.PositiveInfinity;
        }
        private List<Node> VisitedNodes { get; }
        private double Infinity { get; set; }
        private int Index { get; set; }
    }

    public class FloydWarshallAlgorithm
    {
        public double[,] MatrixOfTheShortesPathes { get; private set; }
        public FloydWarshallAlgorithm(Graph Graph)
        {
            MatrixOfTheShortesPathes = new double[Graph.NumberOfNodes, Graph.NumberOfNodes];

            InnerGraph = (Graph)Graph.Clone();

            Infinity = double.PositiveInfinity;

            for (int FirstIndex = 0; FirstIndex < InnerGraph.NumberOfNodes; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < InnerGraph.NumberOfNodes; SecondIndex++)
                {
                    if (InnerGraph.SetOfNodes[FirstIndex].Connections[SecondIndex].Equals(true))
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
                        MatrixOfTheShortesPathes[SecondIndex, ThirdIndex] = Math.Min(MatrixOfTheShortesPathes[SecondIndex, ThirdIndex], MatrixOfTheShortesPathes[SecondIndex, FirstIndex] + MatrixOfTheShortesPathes[FirstIndex, ThirdIndex]);
                    }
                }
            }

            return true;
        }
        private double Infinity { get; set; }
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
            MatrixOfTheShortesPathes = new double[Graph.NumberOfNodes + 1, Graph.NumberOfNodes + 1];

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

    public class FordFulkersonAlgorithm
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

                    CurrentStreamCapacity = Math.Min(CurrentStreamCapacity, InnerGraph.FindEdge(CurrentNode, InnerGraph.SetOfNodes[Index]).Weight);
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

                    CurrentStreamCapacity = Math.Min(CurrentStreamCapacity, InnerGraph.FindEdge(CurrentNode, InnerGraph.SetOfNodes[Index]).Weight);
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

            Infinity = double.PositiveInfinity;

            CapacityOfFlow = 0.0;
        }
        private static Graph InnerGraph { get; set; }
        private Node[] FoundPath { get; set; }
        private double Infinity { get; set; }
    }

    public class AStarAlgorithm
    {
        public bool FindTheShortesPath(Node Start, Node Target)
        {
            OpenSet.Add(Start);

            Node CurrentNode = new Node();

            for (; OpenSet.Count > 0;)
            {
                CurrentNode = FindMinWeightNode();

                OpenSet.Remove(CurrentNode);

                ClosedSet.Add(CurrentNode);

                if(CurrentNode.Equals(Target))
                {
                    return true;
                }
            }

            return true;
        }
        private static Graph InnerGraph { get; set; }
        private List<Node> ClosedSet { get; set; }
        private List<Node> OpenSet { get; set; }
        public AStarAlgorithm(Graph Graph)
        {
            InnerGraph = (Graph)Graph.Clone();

            ClosedSet = new List<Node>();

            OpenSet = new List<Node>();
        }
        private Node FindMinWeightNode()
        {
            Node TemporaryNode = new Node();

            for(int Index = 0; Index < OpenSet.Count; Index++)
            {
                if(TemporaryNode.Weight > OpenSet[Index].Weight)
                {
                    TemporaryNode = OpenSet[Index];
                }
            }

            return TemporaryNode;
        }
    }

    class Program
    {
        private static Graph Graph = new Graph();

        private static Graph NewGraph = new Graph();

        static void Main(string[] args)
        {
            for (int Index = 0; Index < 7; Index++)
            {
                Graph.AddNode(new Node(7, Convert.ToString(Index)));
            }

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

            Console.WriteLine("[Current Graph]\n");

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

            DepthFirstSearchAlgorithm DFS = new DepthFirstSearchAlgorithm(Graph);

            Console.WriteLine($"[1 : Depth First Search]\n\n[All nodes were visited : {DFS.DepthFirstSearch(Graph.SetOfNodes[0])}]\n");

            Console.WriteLine(string.Join(" => ", DFS.ResultOfSearching.Select(Edge => Edge.Name)));

            Console.WriteLine();

            BreadthFirstSearchAlgorithm BFS = new BreadthFirstSearchAlgorithm(Graph);

            Console.WriteLine($"[2 : Breadth First Search]\n\n[All nodes were visited : {BFS.BreadthFirstSearch(Graph.SetOfNodes[0])}]\n");

            Console.WriteLine(string.Join(" => ", BFS.ResultOfSearching.Select(Edge => Edge.Name)));

            Console.WriteLine();

            KruscalAlgorithm KruskalTreeSearch = new KruscalAlgorithm(Graph);

            Console.WriteLine($"[3 : Kruskal Tree Search]\n\n[All nodes were visited : {KruskalTreeSearch.FindMinimumSpanningTree()}]\n");

            Console.WriteLine(string.Join(" => ", KruskalTreeSearch.MinimumSpanningTree.Select(Edges => Edges.Name)));

            Console.WriteLine();

            PrimAlgorithm PrimTreeSearch = new PrimAlgorithm(Graph);

            Console.WriteLine($"[4 : Prim Tree Search]\n\n[All nodes were visited : {PrimTreeSearch.FindMinimumSpanningTree()}]\n");

            Console.WriteLine(string.Join(" => ", PrimTreeSearch.MinimumSpanningTree.Select(Edges => Edges.Name)));

            Console.WriteLine();

            Console.WriteLine("[5 : Bellman-Ford Path Search]\n");

            BellmanFordAlgorithm BellmanFordPathSearch = new BellmanFordAlgorithm(Graph);

            Console.WriteLine($"[Doesn't contain negative cycle : {BellmanFordPathSearch.FindTheShortestPathes(Graph.SetOfNodes[0])}]\n");

            Console.Write($"{Graph.SetOfNodes[0].Name} => ");

            Console.WriteLine(string.Join($"{Graph.SetOfNodes[0].Name} => ", BellmanFordPathSearch.ShortestPathes.Select(Nodes => Nodes.Name + $" = {Nodes.Weight}\n")));

            Console.WriteLine("[6 : Dijkstra Path Search]\n");

            DijkstraAlgorithm DijkstraPathSearch = new DijkstraAlgorithm(Graph);

            Console.WriteLine($"[Doesn't contain negative cycle : {DijkstraPathSearch.FindTheShortestPathes(Graph.SetOfNodes[0])}]\n");

            Console.Write($"{Graph.SetOfNodes[0].Name} => ");

            Console.WriteLine(string.Join($"{Graph.SetOfNodes[0].Name} => ", DijkstraPathSearch.ShortestPathes.Select(Nodes => Nodes.Name + $" = {Nodes.Weight}\n")));

            FloydWarshallAlgorithm FloydWarshallPathSearch = new FloydWarshallAlgorithm(Graph);

            FloydWarshallPathSearch.FindAllTheShortestPathes();

            Console.WriteLine("[7 : Floyd-Warshall Path Search]\n");

            Console.Write("\t");

            for (int Index = 0; Index < Graph.SetOfNodes.Count; Index++)
            {
                Console.Write($"{Graph.SetOfNodes[Index].Name}\t");
            }

            Console.WriteLine("\n");

            for (int FirstIndex = 0; FirstIndex < Graph.SetOfNodes.Count; FirstIndex++)
            {
                Console.Write($"{Graph.SetOfNodes[FirstIndex].Name}\t");

                for (int SecondIndex = 0; SecondIndex < Graph.SetOfNodes.Count; SecondIndex++)
                {
                    Console.Write($"{FloydWarshallPathSearch.MatrixOfTheShortesPathes[FirstIndex, SecondIndex]} \t");
                }

                Console.WriteLine("\n");
            }

            Console.WriteLine();

            JohnsonAlgorithm JohnsonPathSearch = new JohnsonAlgorithm(Graph);

            JohnsonPathSearch.FindAllTheShortestPathes();

            Console.WriteLine("[8 : Johnson Path Search]\n");

            Console.Write("\t");

            for (int Index = 0; Index < Graph.SetOfNodes.Count; Index++)
            {
                Console.Write($"{Graph.SetOfNodes[Index].Name}\t");
            }

            Console.WriteLine("\n");

            for (int FirstIndex = 0; FirstIndex < Graph.SetOfNodes.Count; FirstIndex++)
            {
                Console.Write($"{Graph.SetOfNodes[FirstIndex].Name}\t");

                for (int SecondIndex = 0; SecondIndex < Graph.SetOfNodes.Count; SecondIndex++)
                {
                    Console.Write($"{JohnsonPathSearch.MatrixOfTheShortesPathes[FirstIndex, SecondIndex]} \t");
                }

                Console.WriteLine("\n");
            }

            Console.WriteLine();

            Console.WriteLine("[9 : Ford-Fulkerson Stream Search]\n");

            for (int Index = 0; Index < 7; Index++)
            {
                NewGraph.AddNode(new Node(7, Convert.ToString(Index)));
            }

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

            FordFulkersonAlgorithm FordFulkersonStreamSearch = new FordFulkersonAlgorithm(NewGraph);

            FordFulkersonStreamSearch.FindMaximumFlowWithDFSRealization(NewGraph.SetOfNodes[0], NewGraph.SetOfNodes[NewGraph.NumberOfNodes - 1]);

            Console.WriteLine($"Maximum stream capacity is : {FordFulkersonStreamSearch.CapacityOfFlow}");

            Console.WriteLine("\n[All algorithms were calculated!]");

            Console.ReadKey();
        }
    }
}