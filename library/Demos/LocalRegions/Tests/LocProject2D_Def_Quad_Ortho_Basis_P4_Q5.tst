<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project2D Deformed Quad Ortho basis P=4 Q=5</description>
    <executable>LocProject</executable>
    <parameters>-s quadrilateral -b Ortho_A Ortho_A -o 6 6 -p 7 7 -c 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.0970116</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-11">0.793473</value>
        </metric>
    </metrics>
</test>
