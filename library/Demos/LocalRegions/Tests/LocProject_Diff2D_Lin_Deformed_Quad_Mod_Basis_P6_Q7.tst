<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>LocProject_Diff2D Lin. Deformed Quad Mod. Basis P=6, Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s quadrilateral -b Modified_A Modified_A -o 6 6 -p 7 7 -c 0.0 0.0 1.0 0.0 1.5 1.5 0.0 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">0.288711</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">2.40284</value>
        </metric>
    </metrics>
</test>


