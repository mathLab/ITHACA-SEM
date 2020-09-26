<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project2D Quad Modified basis P=4 Q=5</description>
    <executable>LocProject</executable>
    <parameters>-s quadrilateral -b Modified_A Modified_A -o 6 6 -p 7 7 -c 0.0 0.0 1.0 0.0 1.0 1.0 0.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">1.24683e-13</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">5.82645e-13</value>
        </metric>
    </metrics>
</test>


