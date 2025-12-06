import javax.swing.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.util.List;
import java.util.Queue;

public class Maze extends JFrame {
    private MazePanel mazePanel;
    private JComboBox<String> algorithmSelector;
    private JButton solveButton, resetButton;
    private JLabel statusLabel;

    public Maze() {
        setTitle("Algorithms Maze Pathfinder");
        setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
        setLayout(new BorderLayout(10, 10));

        // Create maze panel
        mazePanel = new MazePanel();
        mazePanel.setParentFrame(this);
        add(mazePanel, BorderLayout.CENTER);

        // Create control panel
        JPanel controlPanel = createControlPanel();
        add(controlPanel, BorderLayout.NORTH);

        // Create legend panel
        JPanel legendPanel = createLegendPanel();
        add(legendPanel, BorderLayout.SOUTH);

        pack();
        setLocationRelativeTo(null);
        setResizable(false);
    }

    private JPanel createControlPanel() {
        JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT, 10, 10));
        panel.setBorder(BorderFactory.createEmptyBorder(10, 10, 10, 10));

        solveButton = new JButton("Solve Maze");
        solveButton.addActionListener(e -> mazePanel.solveMaze());

        JComboBox<String> mazeGeneratorSelector = new JComboBox<>(new String[]{"Kruskal", "Prim"});
        
        resetButton = new JButton("Generate New Maze");
        resetButton.addActionListener(e -> mazePanel.generateNewMaze());

        algorithmSelector = new JComboBox<>(new String[]{"BFS Algorithm", "DFS Algorithm", "Dijkstra Algorithm", "A* Algorithm"});
        algorithmSelector.addActionListener(e -> mazePanel.setAlgorithm(algorithmSelector.getSelectedIndex()));

        statusLabel = new JLabel("Ready");
        statusLabel.setFont(new Font("Arial", Font.BOLD, 12));

        panel.add(new JLabel("Algorithm:"));
        panel.add(algorithmSelector);
        panel.add(solveButton);
        panel.add(mazeGeneratorSelector);
        panel.add(resetButton);
        panel.add(new JLabel("  |  Status:"));
        panel.add(statusLabel);

        return panel;
    }

    private JPanel createLegendPanel() {
        JPanel panel = new JPanel(new FlowLayout(FlowLayout.LEFT, 15, 5));
        panel.setBorder(BorderFactory.createTitledBorder("Legend"));

        addLegendItem(panel, "Start", new Color(76, 175, 80));
        addLegendItem(panel, "Exit", new Color(244, 67, 54));
        addLegendItem(panel, "Normal (0)", Color.WHITE);
        addLegendItem(panel, "Grass (1)", new Color(144, 238, 144));
        addLegendItem(panel, "Mud (5)", new Color(139, 115, 85));
        addLegendItem(panel, "Water (10)", new Color(65, 105, 225));
        addLegendItem(panel, "Visited", new Color(225, 190, 231));
        addLegendItem(panel, "Current", new Color(255, 152, 0));
        addLegendItem(panel, "Path", new Color(255, 215, 0));

        return panel;
    }

    private void addLegendItem(JPanel panel, String label, Color color) {
        JPanel item = new JPanel(new FlowLayout(FlowLayout.LEFT, 5, 0));
        JPanel colorBox = new JPanel();
        colorBox.setPreferredSize(new Dimension(20, 20));
        colorBox.setBackground(color);
        colorBox.setBorder(BorderFactory.createLineBorder(Color.BLACK));
        item.add(colorBox);
        item.add(new JLabel(label));
        panel.add(item);
    }

    public void setStatus(String status) {
        statusLabel.setText(status);
    }

    public static void main(String[] args) {
        SwingUtilities.invokeLater(() -> {
            Maze frame = new Maze();
            frame.setVisible(true);
        });
    }
}

// Cell class representing each maze cell
class Cell {
    int x, y, id;
    TerrainType terrain;
    boolean topWall, rightWall, bottomWall, leftWall;

    public Cell(int x, int y, int id) {
        this.x = x;
        this.y = y;
        this.id = id;
        this.topWall = this.rightWall = this.bottomWall = this.leftWall = true;
        
        // Random terrain assignment
        double rand = Math.random();
        if (rand < 0.2) terrain = TerrainType.GRASS;
        else if (rand < 0.35) terrain = TerrainType.MUD;
        else if (rand < 0.45) terrain = TerrainType.WATER;
        else terrain = TerrainType.NORMAL;
    }
}

// Terrain types with weights
enum TerrainType {
    NORMAL(0, new Color(255, 255, 255)),
    GRASS(1, new Color(144, 238, 144)),
    MUD(5, new Color(139, 115, 85)),
    WATER(10, new Color(65, 105, 225));

    final int weight;
    final Color color;

    TerrainType(int weight, Color color) {
        this.weight = weight;
        this.color = color;
    }
}

// Edge for Kruskal's algorithm
class Edge {
    int from, to;
    String direction;

    public Edge(int from, int to, String direction) {
        this.from = from;
        this.to = to;
        this.direction = direction;
    }
}

// Disjoint Set Union-Find for Kruskal's algorithm
class DisjointSet {
    private int[] parent, rank;

    public DisjointSet(int size) {
        parent = new int[size];
        rank = new int[size];
        for (int i = 0; i < size; i++) {
            parent[i] = i;
            rank[i] = 0;
        }
    }

    public int find(int x) {
        if (parent[x] != x) {
            parent[x] = find(parent[x]);
        }
        return parent[x];
    }

    public boolean union(int x, int y) {
        int rootX = find(x);
        int rootY = find(y);

        if (rootX == rootY) return false;

        if (rank[rootX] < rank[rootY]) {
            parent[rootX] = rootY;
        } else if (rank[rootX] > rank[rootY]) {
            parent[rootY] = rootX;
        } else {
            parent[rootY] = rootX;
            rank[rootX]++;
        }
        return true;
    }
}

// Node for pathfinding algorithms
class PathNode implements Comparable<PathNode> {
    Cell cell;
    double gScore, fScore;

    public PathNode(Cell cell, double gScore, double fScore) {
        this.cell = cell;
        this.gScore = gScore;
        this.fScore = fScore;
    }

    @Override
    public int compareTo(PathNode other) {
        return Double.compare(this.fScore, other.fScore);
    }
}

// Main Maze Panel
class MazePanel extends JPanel {
    private static final int CELL_SIZE = 35;
    private static final int MAZE_WIDTH = 20;
    private static final int MAZE_HEIGHT = 15;
    private static final int ANIMATION_DELAY = 30;

    private Cell[][] maze;
    private Set<String> visitedCells;
    private List<Cell> finalPath;
    private Cell currentCell;
    private int algorithm; // 0 = Dijkstra, 1 = A*
    private boolean isAnimating;
    private Maze parentFrame;

    public MazePanel() {
        setPreferredSize(new Dimension(MAZE_WIDTH * CELL_SIZE + 1, MAZE_HEIGHT * CELL_SIZE + 1));
        setBackground(Color.WHITE);
        visitedCells = new HashSet<>();
        finalPath = new ArrayList<>();
        algorithm = 0;
        generateMazeKruskal();
    }
    
    public void setParentFrame(Maze frame) {
        this.parentFrame = frame;
    }

    public void setAlgorithm(int algo) {
        this.algorithm = algo;
    }

    public void generateNewMaze() {
        if (isAnimating) return;
        visitedCells.clear();
        finalPath.clear();
        currentCell = null;
        
        if (algorithm == 0) {
        generateMazeKruskal();
        } else if (algorithm == 1) {
            generateMazePrim();
        } else {
            generateMazeKruskal();
        }

        repaint();
        if (parentFrame != null) {
            parentFrame.setStatus("New maze generated");
        }
    }

    private void generateMazeKruskal() {
        maze = new Cell[MAZE_HEIGHT][MAZE_WIDTH];
        List<Edge> edges = new ArrayList<>();

        // Initialize cells
        for (int y = 0; y < MAZE_HEIGHT; y++) {
            for (int x = 0; x < MAZE_WIDTH; x++) {
                maze[y][x] = new Cell(x, y, y * MAZE_WIDTH + x);
            }
        }

        // Create all possible edges
        for (int y = 0; y < MAZE_HEIGHT; y++) {
            for (int x = 0; x < MAZE_WIDTH; x++) {
                int id = y * MAZE_WIDTH + x;
                if (x < MAZE_WIDTH - 1) {
                    edges.add(new Edge(id, id + 1, "right"));
                }
                if (y < MAZE_HEIGHT - 1) {
                    edges.add(new Edge(id, id + MAZE_WIDTH, "bottom"));
                }
            }
        }

        // Kruskal's algorithm - shuffle edges
        Collections.shuffle(edges);
        DisjointSet dsu = new DisjointSet(MAZE_WIDTH * MAZE_HEIGHT);

        for (Edge edge : edges) {
            if (dsu.union(edge.from, edge.to)) {
                removeWall(edge);
            }
        }

        // Add extra passages for multiple paths (15% more)
        int extraPassages = (int) (edges.size() * 0.15);
        Collections.shuffle(edges);
        for (int i = 0; i < extraPassages && i < edges.size(); i++) {
            removeWall(edges.get(i));
        }
    }

    private void generateMazePrim() {
        
    maze = new Cell[MAZE_HEIGHT][MAZE_WIDTH];

    // Initialize all cells with walls
    for (int y = 0; y < MAZE_HEIGHT; y++) {
        for (int x = 0; x < MAZE_WIDTH; x++) {
            maze[y][x] = new Cell(x, y, y * MAZE_WIDTH + x);
        }
    }

    boolean[][] visited = new boolean[MAZE_HEIGHT][MAZE_WIDTH];
    List<Edge> frontier = new ArrayList<>();

    // Start at random cell
    int startX = (int)(Math.random() * MAZE_WIDTH);
    int startY = (int)(Math.random() * MAZE_HEIGHT);
    visited[startY][startX] = true;

    // Add walls from start cell
    addFrontierWalls(startX, startY, frontier);

    Random rand = new Random();

        while (!frontier.isEmpty()) {
            // Random frontier wall
            Edge edge = frontier.remove(rand.nextInt(frontier.size()));

            int fromY = edge.from / MAZE_WIDTH;
            int fromX = edge.from % MAZE_WIDTH;
            int toY   = edge.to   / MAZE_WIDTH;
            int toX   = edge.to   % MAZE_WIDTH;

            // Only carve if destination not visited
            if (!visited[toY][toX]) {
                visited[toY][toX] = true;
                removeWall(edge);

                // Add frontier walls of newly visited cell
                addFrontierWalls(toX, toY, frontier);
            }
        }
    }

private void addFrontierWalls(int x, int y, List<Edge> frontier) {
    int id = y * MAZE_WIDTH + x;

    if (x > 0)
        frontier.add(new Edge(id, id - 1, "left"));
    if (x < MAZE_WIDTH - 1)
        frontier.add(new Edge(id, id + 1, "right"));
    if (y > 0)
        frontier.add(new Edge(id, id - MAZE_WIDTH, "top"));
    if (y < MAZE_HEIGHT - 1)
        frontier.add(new Edge(id, id + MAZE_WIDTH, "bottom"));
    }

    private void removeWall(Edge edge) {
        int fromY = edge.from / MAZE_WIDTH;
        int fromX = edge.from % MAZE_WIDTH;
        int toY = edge.to / MAZE_WIDTH;
        int toX = edge.to % MAZE_WIDTH;

        Cell fromCell = maze[fromY][fromX];
        Cell toCell = maze[toY][toX];

        if (edge.direction.equals("right")) {
            fromCell.rightWall = false;
            toCell.leftWall = false;
        } else {
            fromCell.bottomWall = false;
            toCell.topWall = false;
        }
    }

    public void solveMaze() {
        if (isAnimating) return;
        
        visitedCells.clear();
        finalPath.clear();
        currentCell = null;
        isAnimating = true;
        
        if (parentFrame != null) {
            String algoName = switch (algorithm) {
                case 0 -> "BFS";
                case 1 -> "DFS";
                case 2 -> "Dijkstra";
                case 3 -> "A*";
                default -> "Unknown Algorithm";
            };
        parentFrame.setStatus("Solving with " + algoName + "...");
        }

        new Thread(() -> {
            Cell start = maze[0][0];
            Cell end = maze[MAZE_HEIGHT - 1][MAZE_WIDTH - 1];

            if (algorithm == 0) {
                dijkstra(start, end);
            } else if (algorithm == 1){
                aStar(start, end);
            }else if (algorithm == 2){
                bfs(start, end);
            }else{
                dfs(start, end);
            }

            currentCell = null;
            isAnimating = false;
            SwingUtilities.invokeLater(() -> {
                repaint();
                if (parentFrame != null) {
                    parentFrame.setStatus("Path found! Length: " + finalPath.size() + " cells");
                }
            });
        }).start();
    }

    private void bfs(Cell start, Cell end) {
        Queue<Cell> queue = new LinkedList<>();
        Map<String, Cell> previous = new HashMap<>();
        Set<String> visited = new HashSet<>();

        queue.add(start);
        visited.add(start.x + "," + start.y);

        while (!queue.isEmpty()) {
            Cell current = queue.poll();
            String key = current.x + "," + current.y;

            visitedCells.add(key);
            currentCell = current;
            animateDelay(); 
            SwingUtilities.invokeLater(this::repaint);

            if (current == end) {
                reconstructPath(previous, end);
                return;
            }

            for (Cell neighbor : getNeighbors(current)) {
                String nk = neighbor.x + "," + neighbor.y;
                if (!visited.contains(nk)) {
                    visited.add(nk);
                    previous.put(nk, current);
                    queue.add(neighbor);
                }
            }
        }
    }

    private void dfs(Cell start, Cell end) {
        Stack<Cell> stack = new Stack<>();
        Map<String, Cell> previous = new HashMap<>();
        Set<String> visited = new HashSet<>();

        stack.push(start);
        visited.add(start.x + "," + start.y);

        while (!stack.isEmpty()) {
            Cell current = stack.pop();
            String key = current.x + "," + current.y;

            visitedCells.add(key);
            currentCell = current;
            animateDelay(); 
            SwingUtilities.invokeLater(this::repaint);

            if (current == end) {
                reconstructPath(previous, end);
                return;
            }

            for (Cell neighbor : getNeighbors(current)) {
                String nk = neighbor.x + "," + neighbor.y;
                if (!visited.contains(nk)) {
                    visited.add(nk);
                    previous.put(nk, current);
                    stack.push(neighbor);
                }
            }
        }
    }

    private void dijkstra(Cell start, Cell end) {
        Map<String, Double> distances = new HashMap<>();
        Map<String, Cell> previous = new HashMap<>();
        PriorityQueue<PathNode> pq = new PriorityQueue<>();
        Set<String> visited = new HashSet<>();

        // Initialize distances
        for (int y = 0; y < MAZE_HEIGHT; y++) {
            for (int x = 0; x < MAZE_WIDTH; x++) {
                distances.put(x + "," + y, Double.POSITIVE_INFINITY);
            }
        }

        String startKey = start.x + "," + start.y;
        distances.put(startKey, 0.0);
        pq.offer(new PathNode(start, 0, 0));

        while (!pq.isEmpty()) {
            PathNode current = pq.poll();
            String currentKey = current.cell.x + "," + current.cell.y;

            if (visited.contains(currentKey)) continue;
            visited.add(currentKey);

            visitedCells.add(currentKey);
            currentCell = current.cell;
            animateDelay();

            if (current.cell == end) {
                reconstructPath(previous, end);
                return;
            }

            List<Cell> neighbors = getNeighbors(current.cell);
            for (Cell neighbor : neighbors) {
                String neighborKey = neighbor.x + "," + neighbor.y;
                if (visited.contains(neighborKey)) continue;

                double newDist = distances.get(currentKey) + neighbor.terrain.weight + 1;
                if (newDist < distances.get(neighborKey)) {
                    distances.put(neighborKey, newDist);
                    previous.put(neighborKey, current.cell);
                    pq.offer(new PathNode(neighbor, newDist, newDist));
                }
            }
        }
    }

    private void aStar(Cell start, Cell end) {
        Map<String, Double> gScore = new HashMap<>();
        Map<String, Cell> previous = new HashMap<>();
        PriorityQueue<PathNode> openSet = new PriorityQueue<>();
        Set<String> visited = new HashSet<>();

        // Initialize g scores
        for (int y = 0; y < MAZE_HEIGHT; y++) {
            for (int x = 0; x < MAZE_WIDTH; x++) {
                gScore.put(x + "," + y, Double.POSITIVE_INFINITY);
            }
        }

        String startKey = start.x + "," + start.y;
        gScore.put(startKey, 0.0);
        openSet.offer(new PathNode(start, 0, heuristic(start, end)));

        while (!openSet.isEmpty()) {
            PathNode current = openSet.poll();
            String currentKey = current.cell.x + "," + current.cell.y;

            if (visited.contains(currentKey)) continue;
            visited.add(currentKey);

            visitedCells.add(currentKey);
            currentCell = current.cell;
            animateDelay();

            if (current.cell == end) {
                reconstructPath(previous, end);
                return;
            }

            List<Cell> neighbors = getNeighbors(current.cell);
            for (Cell neighbor : neighbors) {
                String neighborKey = neighbor.x + "," + neighbor.y;
                if (visited.contains(neighborKey)) continue;

                double tentativeG = gScore.get(currentKey) + neighbor.terrain.weight + 1;
                if (tentativeG < gScore.get(neighborKey)) {
                    previous.put(neighborKey, current.cell);
                    gScore.put(neighborKey, tentativeG);
                    double fScore = tentativeG + heuristic(neighbor, end);
                    openSet.offer(new PathNode(neighbor, tentativeG, fScore));
                }
            }
        }
    }

    private double heuristic(Cell a, Cell b) {
        return Math.abs(a.x - b.x) + Math.abs(a.y - b.y);
    }

    private List<Cell> getNeighbors(Cell cell) {
        List<Cell> neighbors = new ArrayList<>();
        
        if (!cell.topWall && cell.y > 0)
            neighbors.add(maze[cell.y - 1][cell.x]);
        if (!cell.rightWall && cell.x < MAZE_WIDTH - 1)
            neighbors.add(maze[cell.y][cell.x + 1]);
        if (!cell.bottomWall && cell.y < MAZE_HEIGHT - 1)
            neighbors.add(maze[cell.y + 1][cell.x]);
        if (!cell.leftWall && cell.x > 0)
            neighbors.add(maze[cell.y][cell.x - 1]);

        return neighbors;
    }

    private void reconstructPath(Map<String, Cell> previous, Cell end) {
        finalPath.clear();
        String currentKey = end.x + "," + end.y;
        
        while (currentKey != null) {
            String[] parts = currentKey.split(",");
            int x = Integer.parseInt(parts[0]);
            int y = Integer.parseInt(parts[1]);
            finalPath.add(0, maze[y][x]);
            
            Cell prev = previous.get(currentKey);
            currentKey = prev != null ? prev.x + "," + prev.y : null;
        }
    }
    private void animateDelay() {
        try {
            Thread.sleep(ANIMATION_DELAY);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        SwingUtilities.invokeLater(() -> repaint());
    }

    @Override
    protected void paintComponent(Graphics g) {
        super.paintComponent(g);
        Graphics2D g2d = (Graphics2D) g;
        g2d.setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);

        if (maze == null) return;

        for (int y = 0; y < MAZE_HEIGHT; y++) {
            for (int x = 0; x < MAZE_WIDTH; x++) {
                Cell cell = maze[y][x];
                int px = x * CELL_SIZE;
                int py = y * CELL_SIZE;

                // Determine cell color
                Color cellColor = cell.terrain.color;
                boolean isStart = (x == 0 && y == 0);
                boolean isEnd = (x == MAZE_WIDTH - 1 && y == MAZE_HEIGHT - 1);
                boolean isInPath = finalPath.contains(cell);
                boolean isVisited = visitedCells.contains(x + "," + y);
                boolean isCurrent = (currentCell == cell);

                if (isStart) cellColor = new Color(76, 175, 80);
                else if (isEnd) cellColor = new Color(244, 67, 54);
                else if (isCurrent) cellColor = new Color(255, 152, 0);
                else if (isInPath) cellColor = new Color(255, 215, 0);
                else if (isVisited) cellColor = new Color(225, 190, 231);

                g2d.setColor(cellColor);
                g2d.fillRect(px, py, CELL_SIZE, CELL_SIZE);

                // Draw walls
                g2d.setColor(Color.BLACK);
                g2d.setStroke(new BasicStroke(3));
                if (cell.topWall) g2d.drawLine(px, py, px + CELL_SIZE, py);
                if (cell.rightWall) g2d.drawLine(px + CELL_SIZE, py, px + CELL_SIZE, py + CELL_SIZE);
                if (cell.bottomWall) g2d.drawLine(px, py + CELL_SIZE, px + CELL_SIZE, py + CELL_SIZE);
                if (cell.leftWall) g2d.drawLine(px, py, px, py + CELL_SIZE);

                // Draw labels
                g2d.setFont(new Font("Arial", Font.BOLD, 14));
                if (isStart) {
                    g2d.setColor(Color.WHITE);
                    g2d.drawString("S", px + CELL_SIZE / 2 - 5, py + CELL_SIZE / 2 + 5);
                } else if (isEnd) {
                    g2d.setColor(Color.WHITE);
                    g2d.drawString("E", px + CELL_SIZE / 2 - 5, py + CELL_SIZE / 2 + 5);
                }
            }
        }
        if (finalPath != null && !finalPath.isEmpty()) {
            g2d.setColor(new Color(0, 120, 215));
            g2d.setStroke(new BasicStroke(3));
            for (int i = 0; i < finalPath.size() - 1; i++) {
                Cell c1 = finalPath.get(i);
                Cell c2 = finalPath.get(i + 1);
                int x1 = c1.x * CELL_SIZE + CELL_SIZE / 2;
                int y1 = c1.y * CELL_SIZE + CELL_SIZE / 2;
                int x2 = c2.x * CELL_SIZE + CELL_SIZE / 2;
                int y2 = c2.y * CELL_SIZE + CELL_SIZE / 2;
                g2d.drawLine(x1, y1, x2, y2);
            }
        }
    }
}