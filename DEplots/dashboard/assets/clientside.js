window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        set_fixedrange: function (autorangeOn, figure) {
            figure = JSON.parse(JSON.stringify(figure));
            for (let key of ["", "3"]) {
                let axis = `yaxis${key}`
                figure["layout"][axis]["fixedrange"] = autorangeOn;
                figure["layout"][axis]["autorange"] = false;

                console.log(axis, figure["layout"][axis]["fixedrange"])
            }
            return figure

        },
        autorange_graph: function (relayout, figure, autorangeOn) {
            console.log(autorangeOn)
            if (!autorangeOn) {

                return window.dash_clientside.no_update

            }
            figure = JSON.parse(JSON.stringify(figure));
            console.log(relayout)

            let [lxmin, lxmax] = figure["layout"]["xaxis"]["range"];

            let arr = figure["data"][0]["x"]
            let xmax = arr[arr.length - 1]
            let xmin = arr[0]
            let step = arr[1] - arr[0]
            let xminidx = Math.max(Math.floor((lxmin - xmin) / step), 0)
            let xmaxidx = Math.min(Math.ceil((lxmax - xmin) / step), xmax)
            console.log(xmax, xmin, xmaxidx, xminidx, xmin)
            console.log(arr)
            for (let key of ["", "3"]) {
                let axisName = `x${key}`
                console.log(axisName, key)
                let maxY = 0;

                figure.data.forEach(trace => {
                    console.log(trace.yaxis)
                    // Ensure the slice does not exceed the bounds of the y array
                    if (trace.xaxis === axisName) {
                        let ySlice = trace.y.slice(xminidx, xmaxidx);

                        // Find the maximum value within the sliced y-values
                        let traceMax = Math.max(...ySlice);

                        if (traceMax > maxY) {
                            maxY = traceMax;
                        }
                    }
                });
                maxY = maxY + maxY * 0.01
                let axis = `yaxis${key}`

                figure["layout"][axis]["range"][1] = maxY;
                figure["layout"][axis]["range"][0] = 0 - maxY * 0.01;
                figure["layout"][axis]["autorange"] = false;
                figure["layout"][axis]["fixedrange"] = true;
                console.log(axis, figure["layout"][axis]["range"])

            }
            figure["layout"]["yaxis2"]["range"][0] = relayout["yaxis2.range[0]"]
            figure["layout"]["yaxis2"]["range"][1] = relayout["yaxis2.range[0]"] + 5


            return figure

        }
    }
});
