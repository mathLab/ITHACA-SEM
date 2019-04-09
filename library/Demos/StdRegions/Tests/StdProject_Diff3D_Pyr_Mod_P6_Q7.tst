<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>StdProject_Diff3D Pyramid Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 6 6 6 -p 7 7 6 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">4.99628e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">2.5608e-11</value>
        </metric>
    </metrics>
</test>
