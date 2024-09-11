window.dash_clientside = Object.assign({}, window.dash_clientside, {
    clientside: {
        // highlight_row: function (selected_cell, cs, ce, old_cell) {
        //      // cs and ce will only reset the table highlighting on zoom or anything
        //     if (dash_clientside.callback_context.triggered_id === "coverage-start" || dash_clientside.callback_context.triggered_id === "coverage-end") {
        //         selected_cell = undefined;
        //     }
        //     const table = document.getElementById('gff-table');
        //
        //     if (selected_cell) {
        //         let row = selected_cell["row"]
        //         const nthRow = table.getElementsByTagName('tr')[row + 2];
        //         console.log(nthRow)
        //         if (nthRow) {
        //             const cells = nthRow.getElementsByTagName('td'); // Get all cells in the row
        //             for (let i = 0; i < cells.length; i++) {
        //                 cells[i].setAttribute('style', 'background-color: rgba(var(--bs-primary-rgb), 0.2) !important;');
        //             }
        //
        //         }
        //
        //     }
        //
        //     if (old_cell) {
        //         let old_row = old_cell["row"]
        //         const oldRow = table.getElementsByTagName('tr')[old_row + 2];
        //         if (!selected_cell || selected_cell["row_id"] != old_cell["row_id"])
        //             if (oldRow) {
        //                 const oldCells = oldRow.getElementsByTagName('td'); // Get all cells in the row
        //                 let color;
        //                 if (old_row % 2 === 0) {
        //                     color = "var(--bs-tertiary-bg)"
        //                 } else {
        //                     color = "var(--bs-secondary-bg)"
        //                 }
        //                 for (let i = 0; i < oldCells.length; i++) {
        //                     oldCells[i].setAttribute('style', `background-color: ${color} !important;`);
        //                 }
        //
        //             }
        //
        //     }
        //
        //
        //     return selected_cell
        //
        //
        // },
        update_trace_color: function(color, trace_colors, trace){
            console.log(trace_colors)
            trace_colors[trace] = color
            return trace_colors
        },
        highlight_displayed: function (indices, start, end, contig, table_data){
            let subdata = indices.map(i => table_data[i]);
            let visibleIndices = [];
            for (let i = 0; i < subdata.length; i++) {
                if (subdata[i]["seqid"] == contig && subdata[i]["start"] <= end && subdata[i]["end"] >= start) {
                    visibleIndices.push(i+2);
                }
            }
            console.log(visibleIndices);
            const table = document.getElementById('gff-table');
            const rows = table.getElementsByTagName('tr');
            for (let i = 2; i < rows.length; i++) {
                let cells = rows[i].getElementsByTagName('td');
                let color;
                if (visibleIndices.includes(i)) {
                    if (i % 2 === 0) {
                        color = "rgba(var(--bs-primary-rgb), 0.2)";

                    } else {
                        color = "rgba(var(--bs-primary-rgb), 0.3)";
                    }
                } else if (i % 2 === 0) {
                    color = "var(--bs-tertiary-bg)"
                } else {
                    color = "var(--bs-secondary-bg)"
                }
                console.log(i, color)
                for (let k = 0; k < cells.length; k++) {
                    cells[k].setAttribute('style', `background-color: ${color} !important;`);
                }

            }





        },
        set_fixedrange: function (autorangeOn, figure) {
            figure = JSON.parse(JSON.stringify(figure));
            for (let key of ["", "3"]) {
                let axis = `yaxis${key}`
                figure["layout"][axis]["fixedrange"] = autorangeOn;
                figure["layout"][axis]["autorange"] = false;

            }
            return figure

        },
        autorange_graph: function (relayout, figure, autorangeOn) {
            console.log(autorangeOn)
            if (!autorangeOn) {

                return window.dash_clientside.no_update

            }
            figure = JSON.parse(JSON.stringify(figure));

            let [lxmin, lxmax] = figure["layout"]["xaxis"]["range"];

            let arr = figure["data"][0]["x"]
            let xmax = arr[arr.length - 1]
            let xmin = arr[0]
            let step = arr[1] - arr[0]
            let xminidx = Math.max(Math.floor((lxmin - xmin) / step), 0)
            let xmaxidx = Math.min(Math.ceil((lxmax - xmin) / step), xmax)

            for (let key of ["", "3"]) {
                let axisName = `x${key}`
                let maxY = 0;

                figure.data.forEach(trace => {
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
                if (maxY < 10) {
                    maxY = 10;
                }
                maxY = maxY + maxY * 0.01
                let axis = `yaxis${key}`

                figure["layout"][axis]["range"][1] = maxY;
                figure["layout"][axis]["range"][0] = 0 - maxY * 0.01;
                figure["layout"][axis]["autorange"] = false;
                figure["layout"][axis]["fixedrange"] = true;

            }
            //figure["layout"]["yaxis2"]["range"][0] = relayout["yaxis2.range[0]"]
            //figure["layout"]["yaxis2"]["range"][1] = relayout["yaxis2.range[0]"] + 5


            return figure

        }
    }
});
