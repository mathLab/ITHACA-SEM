<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject_Diff2D Quadrilateral Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s quadrilateral -b Modified_A Modified_A -o 6 6 -p 7 7 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-11">1.2852e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">6.05382e-12</value>
        </metric>
    </metrics>
</test>


