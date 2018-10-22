from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import register_cmap

# defining cyan and yellow colormaps to be used in plotting
cyan = {'red':   ((0.0, 1.00, 1.0),
                  (1.0, 0.28, 1.0)),

        'green': ((0.0, 1.00, 1.0),
                  (1.0, 0.81, 1.0)),

        'blue':  ((0.0, 1.00, 1.0),
                  (1.0, 0.80, 1.0))
        }

yellow = {'red':   ((0.0, 1.00, 1.0),
                  (1.0, 0.99, 1.0)),

        'green': ((0.0, 1.00, 1.0),
                  (1.0, 0.84, 1.0)),

        'blue':  ((0.0, 1.00, 1.0),
                  (1.0, 0.00, 1.0))
        }

cyan = LinearSegmentedColormap('cyan', cyan)
register_cmap(cmap=cyan)

yellow = LinearSegmentedColormap('yellow', yellow)
register_cmap(cmap=yellow)
