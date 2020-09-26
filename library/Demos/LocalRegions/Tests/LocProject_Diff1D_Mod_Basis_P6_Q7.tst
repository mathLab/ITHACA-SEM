<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project1D and Diff Segment Modified basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s segment -b Modified_A -o 6 -p 7 -c 0.0 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">2.4869e-14</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">7.05588e-15</value>
        </metric>
    </metrics>
</test>


