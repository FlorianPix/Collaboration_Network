import plotly.graph_objects as go


def vis(coords, locations, connections, title_text='Co-occurring cites in publication affiliations'):
    fig = go.Figure()

    longitudes = []
    latitudes = []
    cities = []
    for location in locations:
        longitudes.append(coords[location].long)
        latitudes.append(coords[location].lat)
        cities.append(location.city)

    fig.add_trace(go.Scattergeo(
        lon=longitudes,
        lat=latitudes,
        hoverinfo='text',
        text=cities,
        mode='markers',
        marker=dict(
            size=2,
            color='rgb(255, 0, 0)',
            line=dict(
                width=3,
                color='rgba(68, 68, 68, 0)'
            )
        )))

    start_longs = []
    end_longs = []
    start_lats = []
    end_lats = []
    num_traces = 0
    weights = []
    for (location_1, location_2, attr) in connections:
        try:
            start_longs.append(coords[location_1].long)
            end_longs.append(coords[location_2].long)
            start_lats.append(coords[location_1].lat)
            end_lats.append(coords[location_2].lat)
            weights.append(attr['weight'])
            num_traces += 1
        except TypeError:
            print(location_1, location_2)
    max_weight = max(weights)

    for i in range(num_traces):
        fig.add_trace(
            go.Scattergeo(
                lon=[start_longs[i], end_longs[i]],
                lat=[start_lats[i], end_lats[i]],
                mode='lines',
                line=dict(width=1, color='red'),
                opacity=min(float(weights[i] / max_weight), 1.0)
            )
        )

    fig.update_layout(
        title_text=title_text,
        showlegend=False,
        geo=dict(
            scope='world',
            projection_type='orthographic',
            showland=True,
            landcolor='rgb(243, 243, 243)',
            countrycolor='rgb(204, 204, 204)',
        ),
    )

    fig.show()
