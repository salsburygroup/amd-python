import plotly.plotly as py
import plotly.graph_objs as go
import pandas
import argparse

parser = argparse.ArgumentParser(
    description='Make plotly heatmap of vmd hydrogen bonds details file', add_help=False
)

# List all possible user input
inputs = parser.add_argument_group('Input arguments')
inputs.add_argument('-h', '--help', action='help')
inputs.add_argument('-vmd',
                    action='store',
                    dest='vmd_details',
                    help='vmd hbond details file',
                    type=str,
                    required=True
                    )
inputs.add_argument('-title',
                    action='store',
                    dest='title',
                    help='Plot titles',
                    type=str,
                    required=True
                    )

UserInput = parser.parse_args()

hbonds = pandas.read_table(UserInput.vmd_details, header=1, skiprows=0, delim_whitespace=True)
hbonds = hbonds.replace('%','',regex=True)
hbonds['occupancy'] = hbonds['occupancy'].apply(pandas.to_numeric)

data = [
    go.Heatmap(
        z = hbonds['occupancy'],
        x = hbonds['donor'],
        y = hbonds['acceptor'],
        colorbar=dict(
            title='% Occupancy',
            titlefont=dict(
                family='Courier New, monospace',
                size=18,
                color='#7f7f7f'
            )
        )
    )
]

layout = go.Layout(
    title=UserInput.title,
    xaxis=dict(
        title='Donor',
        titlefont=dict(
            family='Courier New, monospace',
            size=18,
            color='#7f7f7f'
        )
    ),
    yaxis=dict(
        title='Acceptor',
        titlefont=dict(
            family='Courier New, monospace',
            size=18,
            color='#7f7f7f'
        )
    ),
)
fig = go.Figure(data=data, layout=layout)
py.iplot(fig, filename=UserInput.title)