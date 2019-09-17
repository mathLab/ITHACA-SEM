<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>StdProject_Diff1D Segment Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s segment -b Modified_A -o 6 -p 7 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.32975e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">3.90799e-14</value>
        </metric>
    </metrics>
</test>


