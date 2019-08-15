from matplotlib.colors import LinearSegmentedColormap
from matplotlib.pyplot import register_cmap

# defining cyan and yellow colormaps to be used in plotting
cyan = {
	'red':(
		(0.0, 1.0, 1.0),
		(1.0, 0.0, 1.0)
	),

	'green':(
		(0.0, 1.00000000, 1.00000000),
		(1.0, 0.69019608, 0.69019608)
	),

	'blue':(
		(0.0, 1.00000000, 1.00000000),
		(1.0, 0.94117647, 0.94117647)
	),

	'alpha':(
		(0.0, 0.0, 0.0),
		(1.0, 0.8, 0.8)
	)
}

yellow = {
	'red':(
		(0.0, 1.0, 1.0),
		(1.0, 1.0, 1.0)
	),

	'green':(
		(0.0, 1.00000000, 1.00000000),
		(1.0, 0.75294118, 0.75294118)
	),

	'blue':(
		(0.0, 1.0, 1.0),
		(1.0, 0.0, 0.0)
	),

	'alpha':(
		(0.0, 0.0, 0.0),
		(1.0, 0.8, 0.8)
	)
}


cyan = LinearSegmentedColormap('cyan', cyan)
register_cmap(cmap=cyan)

yellow = LinearSegmentedColormap('yellow', yellow)
register_cmap(cmap=yellow)


class colors:

	reset='\033[0m'
	bold='\033[01m'
	disable='\033[02m'
	underline='\033[04m'
	reverse='\033[07m'
	strikethrough='\033[09m'
	invisible='\033[08m'

	black='\033[30m'
	red='\033[31m'
	green='\033[32m'
	orange='\033[33m'
	blue='\033[34m'
	purple='\033[35m'
	cyan='\033[36m'

	lightgrey='\033[37m'
	darkgrey='\033[90m'
	lightred='\033[91m'
	lightgreen='\033[92m'
	yellow='\033[93m'
	lightblue='\033[94m'
	pink='\033[95m'
	lightcyan='\033[96m'
