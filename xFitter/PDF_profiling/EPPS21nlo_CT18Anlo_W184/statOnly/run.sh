#First, perform profiling for all study cases.
for dir in */; do
    cd ${dir}
    sh run.sh
    cd ..
done

#Some plots depend on results in many dirs, depending on the plot.sh scripts
for dir in */; do
    cd ${dir}
    sh plot.sh
    cd ..
done
