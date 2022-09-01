function  timeAve(halfWindow::Int64 , timeSeries::Array)

    timeAveraged = timeSeries 
    for i = 1 + halfWindow : length(timeSeries) - halfWindow
        timeAveraged[i] = sum(timeSeries[i - halfWindow : i + halfWindow])/(2*halfWindow + 1)
    end
    return timeAveraged
end