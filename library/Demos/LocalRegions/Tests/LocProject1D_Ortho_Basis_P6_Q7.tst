<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project1D Segment Orthogonal basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s Segment -b Ortho_A -o 6 -p 7 -c 0.0 1.0</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">4.26337e-16</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">8.88178e-16</value>
        </metric>
    </metrics>
</test>


