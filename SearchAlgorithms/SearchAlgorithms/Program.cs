using System.Collections.Generic;
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
        public Node(int NumberOfNodes, string NameOfNode)
        {
            Inheritors = new bool[NumberOfNodes + 1];

            this.NameOfNode = $"[{NameOfNode}]";

            WeightOfNode = 0.0;

            IndexOfNode = 0;
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

            WeightOfNode = 0.0;

            IndexOfNode = 0;

            Inheritors = null;
        }
    }

    public class Edge : ICloneable
    {
        public Edge(double WeightOfEdge, params Node[] Nodes)
        {
            EndsOfEdge = new List<Node>();

            for (int Index = 0; Index < Nodes.Length; Index++)
            {
                EndsOfEdge.Add(Nodes[Index]);
            }

            this.WeightOfEdge = WeightOfEdge;

            EndsOfEdge[0].Inheritors[Nodes[1].IndexOfNode] = true;

            EndsOfEdge[1].Inheritors[Nodes[0].IndexOfNode] = true;

            NameOfEdge = $"{Nodes[0].NameOfNode}{Nodes[1].NameOfNode}";
        }
        public List<Node> EndsOfEdge { get; private set; }
        public string NameOfEdge { get; private set; }
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

                NameOfEdge = this.NameOfEdge,

                IndexOfEdge = this.IndexOfEdge,

                EndsOfEdge = EndsOfEdge
            };

            return Edge;
        }
        public Edge()
        {
            EndsOfEdge = new List<Node>();

            NameOfEdge = string.Empty;

            WeightOfEdge = 0.0;

            IndexOfEdge = 0;
        }
    }

    public class Graph : ICloneable
    {
        public Edge GetEdge(Node TheBeginningOfEdge, Node TheEndOfEdge)
        {
            for (int Index = 0; Index < NumberOfEdges; Index++)
            {
                if (SetOfEdges[Index][0].Equals(TheBeginningOfEdge) && SetOfEdges[Index][1].Equals(TheEndOfEdge) || SetOfEdges[Index][0].Equals(TheEndOfEdge) && SetOfEdges[Index][1].Equals(TheBeginningOfEdge))
                {
                    return SetOfEdges[Index];
                }
            }

            return new Edge { WeightOfEdge = double.PositiveInfinity };
        }
        public bool DoesntContainNegativeWeightCycle(List<Node> Nodes)
        {
            for (int Index = 0; Index < NumberOfEdges; Index++)
            {
                if (Nodes[SetOfEdges[Index][1].IndexOfNode].WeightOfNode > Nodes[SetOfEdges[Index][0].IndexOfNode].WeightOfNode + SetOfEdges[Index].WeightOfEdge)
                {
                    return false;
                }
            }
            return true;
        }
        public bool DoesntContainNegativeWeightCycle(Node[] Nodes)
        {
            for (int Index = 0; Index < NumberOfEdges; Index++)
            {
                if (Nodes[SetOfEdges[Index][1].IndexOfNode].WeightOfNode > Nodes[SetOfEdges[Index][0].IndexOfNode].WeightOfNode + SetOfEdges[Index].WeightOfEdge)
                {
                    return false;
                }
            }
            return true;
        }
        public int NumberOfNodes { get { return SetOfNodes.Count; } }
        public int NumberOfEdges { get { return SetOfEdges.Count; } }
        public List<Node> SetOfNodes { get; private set; }
        public List<Edge> SetOfEdges { get; private set; }
        public void AddOneWayEdge(Edge Edge)
        {
            SetOfEdges.Add(Edge);

            SetOfEdges[SetOfEdges.Count - 1].IndexOfEdge = SetOfEdges.Count - 1;
        }
        public void AddTwoWayEdge(Edge Edge)
        {
            SetOfEdges.Add(Edge);

            SetOfEdges[SetOfEdges.Count - 1].IndexOfEdge = SetOfEdges.Count - 1;

            SetOfEdges.Add(new Edge(Edge.WeightOfEdge, Edge[1], Edge[0]));

            SetOfEdges[SetOfEdges.Count - 1].IndexOfEdge = SetOfEdges.Count - 1;
        }
        public void AddNode(Node Node)
        {
            SetOfNodes.Add(Node);

            SetOfNodes[SetOfNodes.Count - 1].IndexOfNode = SetOfNodes.Count - 1;
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
        public bool DFS(Node TheBeginningOfSearching, Node TheEndOfSearching)
        {
            VisitedNodes.Add(TheBeginningOfSearching);

            if(VisitedNodes.Contains(TheEndOfSearching))
            {
                return true;
            }

            for (int Index = 0; Index < AuxiliaryGraph.NumberOfNodes; Index++)
            {
                if (AuxiliaryGraph.SetOfNodes[TheBeginningOfSearching.IndexOfNode].Inheritors[Index].Equals(true) && !VisitedNodes.Contains(AuxiliaryGraph.SetOfNodes[Index]))
                {
                    if (DFS(AuxiliaryGraph.SetOfNodes[Index], TheEndOfSearching))
                    {
                        return true;
                    }
                }
            }

            return false;
        }
        public DepthFirstSearchAlgorithm(Graph Graph)
        {
            AuxiliaryGraph = (Graph)Graph.Clone();

            VisitedNodes = new List<Node>();
        }
        public bool DFS(Node TheBeginningOfSearch)
        {
            VisitedNodes.Add(TheBeginningOfSearch);

            for (int Index = 0; Index < AuxiliaryGraph.NumberOfNodes; Index++)
            {
                if (AuxiliaryGraph.SetOfNodes[TheBeginningOfSearch.IndexOfNode].Inheritors[Index].Equals(true) && !VisitedNodes.Contains(AuxiliaryGraph.SetOfNodes[Index]))
                {
                    DFS(AuxiliaryGraph.SetOfNodes[Index]);
                }
            }

            return VisitedNodes.Count.Equals(AuxiliaryGraph.NumberOfNodes);
        }
        private static Graph AuxiliaryGraph { get; set; }
        public List<Node> VisitedNodes { get; }
    }

    public class BreadthFirstSearchAlgorithm
    {
        public bool BFS(Node TheBeginningOfSearching, Node TheEndOfSearching)
        {
            QueuedNodes.Enqueue(TheBeginningOfSearching);

            VisitedNodes.Add(TheBeginningOfSearching);

            for (; QueuedNodes.Count > 0;)
            {
                TheBeginningOfSearching = QueuedNodes.Dequeue();

                for (int Index = 0; Index < AuxiliaryGraph.NumberOfNodes; Index++)
                {
                    if (AuxiliaryGraph.SetOfNodes[TheBeginningOfSearching.IndexOfNode].Inheritors[Index].Equals(true) && !VisitedNodes.Contains(AuxiliaryGraph.SetOfNodes[Index]))
                    {
                        if (VisitedNodes.Contains(TheEndOfSearching))
                        {
                            return true;
                        }

                        QueuedNodes.Enqueue(AuxiliaryGraph.SetOfNodes[Index]);

                        VisitedNodes.Add(AuxiliaryGraph.SetOfNodes[Index]);
                    }
                }
            }

            return false;
        }
        public BreadthFirstSearchAlgorithm(Graph Graph)
        {
            QueuedNodes = new Queue<Node>();

            AuxiliaryGraph = (Graph)Graph.Clone();

            VisitedNodes = new List<Node>();
        }
        public bool BFS(Node TheBeginningOfSearch)
        {
            QueuedNodes.Enqueue(TheBeginningOfSearch);

            VisitedNodes.Add(TheBeginningOfSearch);

            for (; QueuedNodes.Count > 0;)
            {
                TheBeginningOfSearch = QueuedNodes.Dequeue();

                for (int Index = 0; Index < AuxiliaryGraph.NumberOfNodes; Index++)
                {
                    if (AuxiliaryGraph.SetOfNodes[TheBeginningOfSearch.IndexOfNode].Inheritors[Index].Equals(true) && !VisitedNodes.Contains(AuxiliaryGraph.SetOfNodes[Index]))
                    {
                        QueuedNodes.Enqueue(AuxiliaryGraph.SetOfNodes[Index]);

                        VisitedNodes.Add(AuxiliaryGraph.SetOfNodes[Index]);
                    }
                }
            }

            return VisitedNodes.Count.Equals(AuxiliaryGraph.NumberOfNodes);
        }
        private static Graph AuxiliaryGraph { get; set; }
        private Queue<Node> QueuedNodes { get; }
        public List<Node> VisitedNodes { get; }
    }

    public class KruscalAlgorithm
    {
        public List<Edge> MinimumSpanningTree { get; }
        private static Graph AuxiliaryGraph { get; set; }
        private Node FindSubTree(Node Node)
        {
            Node RootOfTree = Node;

            for (; !RootOfTree.Equals(AuxiliaryGraph.SetOfNodes[RootOfTree.IndexOfNode]);)
            {
                RootOfTree = AuxiliaryGraph.SetOfNodes[RootOfTree.IndexOfNode];
            }

            for (; !Node.Equals(RootOfTree);)
            {
                Node OldNode = AuxiliaryGraph.SetOfNodes[Node.IndexOfNode];

                AuxiliaryGraph.SetOfNodes[Node.IndexOfNode] = RootOfTree;

                Node = OldNode;
            }

            return RootOfTree;
        }
        public KruscalAlgorithm(Graph Graph)
        {
            MinimumSpanningTree = new List<Edge>();

            AuxiliaryGraph = (Graph)Graph.Clone();
        }
        public bool KruscalTreeSearch()
        {
            SortGraph();

            for (int Index = 0; Index < AuxiliaryGraph.NumberOfEdges; Index++)
            {
                Node StartNode = FindSubTree(AuxiliaryGraph.SetOfEdges[Index][0]);

                Node EndNode = FindSubTree(AuxiliaryGraph.SetOfEdges[Index][1]);

                if (!StartNode.Equals(EndNode))
                {
                    MinimumSpanningTree.Add(AuxiliaryGraph.SetOfEdges[Index]);

                    AuxiliaryGraph.SetOfNodes[EndNode.IndexOfNode] = StartNode;
                }
            }

            return true;
        }
        private void SortGraph()
        {
            Edge TemporaryEdge = new Edge();

            for (int FirstIndex = 0; FirstIndex < AuxiliaryGraph.NumberOfEdges; FirstIndex++)
            {
                for (int SecondIndex = FirstIndex + 1; SecondIndex < AuxiliaryGraph.NumberOfEdges; SecondIndex++)
                {
                    if (AuxiliaryGraph.SetOfEdges[FirstIndex].WeightOfEdge > AuxiliaryGraph.SetOfEdges[SecondIndex].WeightOfEdge)
                    {
                        TemporaryEdge = AuxiliaryGraph.SetOfEdges[FirstIndex];

                        AuxiliaryGraph.SetOfEdges[FirstIndex] = AuxiliaryGraph.SetOfEdges[SecondIndex];

                        AuxiliaryGraph.SetOfEdges[SecondIndex] = TemporaryEdge;
                    }
                }
            }
        }
    }

    public class PrimAlgorithm
    {
        public List<Node> VisitedNodes { get; private set; }
        public List<Edge> MinimumSpanningTree { get; }
        private static Graph AuxiliaryGraph { get; set; }
        public PrimAlgorithm(Graph Graph)
        {
            MinimumSpanningTree = new List<Edge>();

            AuxiliaryGraph = (Graph)Graph.Clone();

            VisitedNodes = new List<Node>();
        }
        public bool PrimTreeSearch()
        {
            VisitedNodes.Add(AuxiliaryGraph.SetOfNodes[0]);

            AuxiliaryGraph.SetOfNodes.RemoveAt(VisitedNodes[0].IndexOfNode);

            for (; AuxiliaryGraph.NumberOfNodes > 0;)
            {
                int IndexOfTheEasiestEdge = -1;

                for (int Index = 0; Index < AuxiliaryGraph.NumberOfEdges; Index++)
                {
                    if (!VisitedNodes.IndexOf(AuxiliaryGraph.SetOfEdges[Index][0]).Equals(-1) && !AuxiliaryGraph.SetOfNodes.IndexOf(AuxiliaryGraph.SetOfEdges[Index][1]).Equals(-1) || !VisitedNodes.IndexOf(AuxiliaryGraph.SetOfEdges[Index][1]).Equals(-1) && !AuxiliaryGraph.SetOfNodes.IndexOf(AuxiliaryGraph.SetOfEdges[Index][0]).Equals(-1))
                    {
                        if (!IndexOfTheEasiestEdge.Equals(-1))
                        {
                            if (AuxiliaryGraph.SetOfEdges[Index].WeightOfEdge < AuxiliaryGraph.SetOfEdges[IndexOfTheEasiestEdge].WeightOfEdge)
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

                if (!VisitedNodes.IndexOf(AuxiliaryGraph.SetOfEdges[IndexOfTheEasiestEdge][0]).Equals(-1))
                {
                    VisitedNodes.Add(AuxiliaryGraph.SetOfEdges[IndexOfTheEasiestEdge][1]);

                    AuxiliaryGraph.SetOfNodes.Remove(AuxiliaryGraph.SetOfEdges[IndexOfTheEasiestEdge][1]);
                }
                else
                {
                    VisitedNodes.Add(AuxiliaryGraph.SetOfEdges[IndexOfTheEasiestEdge][0]);

                    AuxiliaryGraph.SetOfNodes.Remove(AuxiliaryGraph.SetOfEdges[IndexOfTheEasiestEdge][0]);
                }

                MinimumSpanningTree.Add(AuxiliaryGraph.SetOfEdges[IndexOfTheEasiestEdge]);

                AuxiliaryGraph.SetOfEdges.RemoveAt(IndexOfTheEasiestEdge);
            }

            return true;
        }
    }

    public class BellmanFordAlgorithm
    {
        public bool BellmanFordPathSearch(Node TheBeginningOfSearch)
        {
            for (int Index = 0; Index < AuxiliaryGraph.NumberOfNodes; Index++)
            {
                Pathes[Index].WeightOfNode = Infinity;
            }

            Pathes[TheBeginningOfSearch.IndexOfNode].WeightOfNode = 0;

            for (int FirstIndex = 1; FirstIndex < AuxiliaryGraph.NumberOfNodes - 1; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < AuxiliaryGraph.NumberOfEdges; SecondIndex++)
                {
                     if (Pathes[AuxiliaryGraph.SetOfEdges[SecondIndex][1].IndexOfNode].WeightOfNode > Pathes[AuxiliaryGraph.SetOfEdges[SecondIndex][0].IndexOfNode].WeightOfNode + AuxiliaryGraph.SetOfEdges[SecondIndex].WeightOfEdge)
                     {
                        Pathes[AuxiliaryGraph.SetOfEdges[SecondIndex][1].IndexOfNode].WeightOfNode = Pathes[AuxiliaryGraph.SetOfEdges[SecondIndex][0].IndexOfNode].WeightOfNode + AuxiliaryGraph.SetOfEdges[SecondIndex].WeightOfEdge;
                     }
                }
            }

            return AuxiliaryGraph.DoesntContainNegativeWeightCycle(Pathes);
        }
        private static Graph AuxiliaryGraph { get; set; }
        public BellmanFordAlgorithm(Graph Graph)
        {
            AuxiliaryGraph = (Graph)Graph.Clone();

            Pathes = AuxiliaryGraph.SetOfNodes;

            Infinity = double.PositiveInfinity;
        }
        public List<Node> Pathes { get; private set; }
        private double Infinity { get; set; }
    }

    public class DijkstraAlgorithm
    {
        public bool DijkstraPathSearch(Node TheBeginningOfSearch)
        {
            for (int Index = 0; Index < AuxiliaryGraph.NumberOfNodes; Index++)
            {
                Pathes[Index].WeightOfNode = Infinity;
            }

            Pathes[TheBeginningOfSearch.IndexOfNode].WeightOfNode = 0;

            for (int FirstIndex = 0; FirstIndex < AuxiliaryGraph.NumberOfNodes - 1; FirstIndex++)
            {
                Node TemporaryNode = new Node { WeightOfNode = Infinity };

                TemporaryNode.WeightOfNode = Infinity;

                for (int SecondIndex = 0; SecondIndex < AuxiliaryGraph.NumberOfNodes; SecondIndex++)
                {
                    if (!VisitedNodes.Contains(Pathes[SecondIndex]) && Pathes[SecondIndex].WeightOfNode <= TemporaryNode.WeightOfNode)
                    {
                        TemporaryNode = Pathes[SecondIndex];

                        Index = Pathes[SecondIndex].IndexOfNode;
                    }
                }

                VisitedNodes.Add(Pathes[Index]);

                for (int SecondIndex = 0; SecondIndex < AuxiliaryGraph.NumberOfNodes; SecondIndex++)
                {
                    if (!VisitedNodes.Contains(Pathes[SecondIndex]) && AuxiliaryGraph.SetOfNodes[Index].Inheritors[SecondIndex].Equals(true) && !Pathes[Index].WeightOfNode.Equals(Infinity) && Pathes[SecondIndex].WeightOfNode > Pathes[Index].WeightOfNode + AuxiliaryGraph.GetEdge(AuxiliaryGraph.SetOfNodes[Index], AuxiliaryGraph.SetOfNodes[SecondIndex]).WeightOfEdge)
                    {
                        Pathes[SecondIndex].WeightOfNode = Pathes[Index].WeightOfNode + AuxiliaryGraph.GetEdge(AuxiliaryGraph.SetOfNodes[Index], AuxiliaryGraph.SetOfNodes[SecondIndex]).WeightOfEdge;
                    }
                }
            }

            return AuxiliaryGraph.DoesntContainNegativeWeightCycle(Pathes);
        }
        private static Graph AuxiliaryGraph { get; set; }
        public List<Node> Pathes { get; private set; }
        private List<Node> VisitedNodes { get; }
        public DijkstraAlgorithm(Graph Graph)
        {
            AuxiliaryGraph = (Graph)Graph.Clone();

            Pathes = AuxiliaryGraph.SetOfNodes;

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
            MatrixOfShortesPathes = new double[Graph.NumberOfNodes, Graph.NumberOfNodes];

            AuxiliaryGraph = (Graph)Graph.Clone();

            Infinity = double.PositiveInfinity;

            for (int FirstIndex = 0; FirstIndex < AuxiliaryGraph.NumberOfNodes; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < AuxiliaryGraph.NumberOfNodes; SecondIndex++)
                {
                    if (AuxiliaryGraph.SetOfNodes[FirstIndex].Inheritors[SecondIndex].Equals(true))
                    {
                        MatrixOfShortesPathes[FirstIndex, SecondIndex] = AuxiliaryGraph.GetEdge(AuxiliaryGraph.SetOfNodes[FirstIndex], AuxiliaryGraph.SetOfNodes[SecondIndex]).WeightOfEdge;
                    }
                    else
                    {
                        MatrixOfShortesPathes[FirstIndex, SecondIndex] = Infinity;
                    }
                }
            }
        }
        private static Graph AuxiliaryGraph { get; set; }
        public void FloydWarshallPathSearch()
        {
            for (int Index = 0; Index < AuxiliaryGraph.NumberOfNodes; Index++)
            {
                MatrixOfShortesPathes[Index, Index] = 0;
            }

            for (int FirstIndex = 0; FirstIndex < AuxiliaryGraph.NumberOfNodes; FirstIndex++)
            {
                for (int SecondIndex = 0; SecondIndex < AuxiliaryGraph.NumberOfNodes; SecondIndex++)
                {
                    for (int ThirdIndex = 0; ThirdIndex < AuxiliaryGraph.NumberOfNodes; ThirdIndex++)
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

    public class JohnsonAlgorithm
    {
        private static BellmanFordAlgorithm BellmanFordPathSearch { get; set; }
        private static DijkstraAlgorithm DijkstraPathSearch { get; set; }
        public double[,] MatrixOfShortesPathes { get; private set; }
        private static Graph CopyOfOriginalGraph { get; set; }
        private static Graph AuxiliaryGraph { get; set; }
        public JohnsonAlgorithm(Graph Graph)
        {
            MatrixOfShortesPathes = new double[Graph.NumberOfNodes + 1, Graph.NumberOfNodes + 1];

            AuxiliaryGraph = (Graph)Graph.Clone();

            CopyOfOriginalGraph = (Graph)Graph.Clone();
        }
        public bool JohnsonPathSearch()
        {
            AuxiliaryGraph.AddNode(new Node(AuxiliaryGraph.NumberOfNodes,"Temporary"));

            for(int Index = 0; Index < AuxiliaryGraph.NumberOfNodes - 1; Index++)
            {
                AuxiliaryGraph.AddTwoWayEdge(new Edge(0,AuxiliaryGraph.SetOfNodes[AuxiliaryGraph.NumberOfNodes - 1],AuxiliaryGraph.SetOfNodes[Index]));
            }

            BellmanFordPathSearch = new BellmanFordAlgorithm(AuxiliaryGraph);

            if(BellmanFordPathSearch.BellmanFordPathSearch(AuxiliaryGraph.SetOfNodes[AuxiliaryGraph.NumberOfNodes - 1]))
            {
                for(int FirstIndex = 0; FirstIndex < CopyOfOriginalGraph.NumberOfNodes; FirstIndex++)
                {
                    DijkstraPathSearch = new DijkstraAlgorithm(CopyOfOriginalGraph);

                    DijkstraPathSearch.DijkstraPathSearch(CopyOfOriginalGraph.SetOfNodes[FirstIndex]);

                    for(int SecondIndex = 0; SecondIndex < CopyOfOriginalGraph.NumberOfNodes; SecondIndex++)
                    {
                        MatrixOfShortesPathes[FirstIndex, SecondIndex] = DijkstraPathSearch.Pathes[SecondIndex].WeightOfNode;
                    }
                }

                return true;
            }
            else
            {
                for(int Index = 0; Index < CopyOfOriginalGraph.NumberOfEdges; Index++)
                {
                    CopyOfOriginalGraph.SetOfEdges[Index].WeightOfEdge = CopyOfOriginalGraph.SetOfEdges[Index].WeightOfEdge + BellmanFordPathSearch.Pathes[CopyOfOriginalGraph.SetOfEdges[Index][0].IndexOfNode].WeightOfNode - BellmanFordPathSearch.Pathes[CopyOfOriginalGraph.SetOfEdges[Index][1].IndexOfNode].WeightOfNode;
                }

                for(int FirstIndex = 0; FirstIndex < CopyOfOriginalGraph.NumberOfNodes; FirstIndex++)
                {
                    DijkstraPathSearch = new DijkstraAlgorithm(CopyOfOriginalGraph);

                    DijkstraPathSearch.DijkstraPathSearch(CopyOfOriginalGraph.SetOfNodes[FirstIndex]);

                    for (int SecondIndex = 0; SecondIndex < CopyOfOriginalGraph.NumberOfNodes; SecondIndex++)
                    {
                        MatrixOfShortesPathes[FirstIndex, SecondIndex] = DijkstraPathSearch.Pathes[SecondIndex].WeightOfNode;
                    }
                }
            }

            return true;
        }
    }

    public class FordFulkersonalgorithm
    { }

    class Program
    {
        private static Graph Graph = new Graph();
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
            Console.WriteLine("           11");

            DepthFirstSearchAlgorithm DFS = new DepthFirstSearchAlgorithm(Graph);

            Console.WriteLine($"[1 : Depth First Search]\n\n[All nodes were visited : {DFS.DFS(Graph.SetOfNodes[0])}]\n");

            Console.WriteLine(string.Join(" => ", DFS.VisitedNodes.Select(Vertex => Vertex.NameOfNode)));

            Console.WriteLine();

            BreadthFirstSearchAlgorithm BFS = new BreadthFirstSearchAlgorithm(Graph);

            Console.WriteLine($"[2 : Breadth First Search]\n\n[All nodes were visited : {BFS.BFS(Graph.SetOfNodes[0])}]\n");

            Console.WriteLine(string.Join(" => ", BFS.VisitedNodes.Select(Vertex => Vertex.NameOfNode)));

            Console.WriteLine();

            KruscalAlgorithm KruskalTreeSearch = new KruscalAlgorithm(Graph);

            Console.WriteLine($"[3 : Kruskal Tree Search]\n\n[All nodes were visited : {KruskalTreeSearch.KruscalTreeSearch()}]\n");

            Console.WriteLine(string.Join(" => ", KruskalTreeSearch.MinimumSpanningTree.Select(Edges => Edges.NameOfEdge)));

            Console.WriteLine();

            PrimAlgorithm PrimTreeSearch = new PrimAlgorithm(Graph);

            Console.WriteLine($"[4 : Prim Tree Search]\n\n[All nodes were visited : {PrimTreeSearch.PrimTreeSearch()}]\n");

            Console.WriteLine(string.Join(" => ", PrimTreeSearch.MinimumSpanningTree.Select(Edges => Edges.NameOfEdge)));

            Console.WriteLine();

            Console.WriteLine("[5 : Bellman-Ford Path Search]\n");

            BellmanFordAlgorithm BellmanFordPathSearch = new BellmanFordAlgorithm(Graph);

            Console.WriteLine($"[Doesn't contain negative cycle : {BellmanFordPathSearch.BellmanFordPathSearch(Graph.SetOfNodes[0])}]\n");

            Console.Write($"{Graph.SetOfNodes[0].NameOfNode} => ");

            Console.WriteLine(string.Join($"{Graph.SetOfNodes[0].NameOfNode} => ", BellmanFordPathSearch.Pathes.Select(Nodes => Nodes.NameOfNode + $" = {Nodes.WeightOfNode}\n")));

            Console.WriteLine("[6 : Dijkstra Path Search]\n");

            DijkstraAlgorithm DijkstraPathSearch = new DijkstraAlgorithm(Graph);

            Console.WriteLine($"[Doesn't contain negative cycle : {DijkstraPathSearch.DijkstraPathSearch(Graph.SetOfNodes[0])}]\n");

            Console.Write($"{Graph.SetOfNodes[0].NameOfNode} => ");

            Console.WriteLine(string.Join($"{Graph.SetOfNodes[0].NameOfNode} => ", DijkstraPathSearch.Pathes.Select(Nodes => Nodes.NameOfNode + $" = {Nodes.WeightOfNode}\n")));

            FloydWarshallAlgorithm FloydWarshallPathSearch = new FloydWarshallAlgorithm(Graph);

            FloydWarshallPathSearch.FloydWarshallPathSearch();

            Console.WriteLine("[7 : Floyd-Warshall Path Search]\n");

            Console.Write("\t");

            for (int Index = 0; Index < Graph.SetOfNodes.Count; Index++)
            {
                Console.Write($"{Graph.SetOfNodes[Index].NameOfNode}\t");
            }

            Console.WriteLine("\n");

            for (int FirstIndex = 0; FirstIndex < Graph.SetOfNodes.Count; FirstIndex++)
            {
                Console.Write($"{Graph.SetOfNodes[FirstIndex].NameOfNode}\t");

                for (int SecondIndex = 0; SecondIndex < Graph.SetOfNodes.Count; SecondIndex++)
                {
                    Console.Write($"{FloydWarshallPathSearch.MatrixOfShortesPathes[FirstIndex, SecondIndex]} \t");
                }

                Console.WriteLine("\n");
            }

            Console.WriteLine();

            JohnsonAlgorithm JohnsonPathSearch = new JohnsonAlgorithm(Graph);

            JohnsonPathSearch.JohnsonPathSearch();

            Console.WriteLine("[8 : Johnson Path Search]\n");

            Console.Write("\t");

            for (int Index = 0; Index < Graph.SetOfNodes.Count; Index++)
            {
                Console.Write($"{Graph.SetOfNodes[Index].NameOfNode}\t");
            }

            Console.WriteLine("\n");

            for (int FirstIndex = 0; FirstIndex < Graph.SetOfNodes.Count; FirstIndex++)
            {
                Console.Write($"{Graph.SetOfNodes[FirstIndex].NameOfNode}\t");

                for (int SecondIndex = 0; SecondIndex < Graph.SetOfNodes.Count; SecondIndex++)
                {
                    Console.Write($"{JohnsonPathSearch.MatrixOfShortesPathes[FirstIndex, SecondIndex]} \t");
                }

                Console.WriteLine("\n");
            }

            Console.ReadKey();
        }
    }
}