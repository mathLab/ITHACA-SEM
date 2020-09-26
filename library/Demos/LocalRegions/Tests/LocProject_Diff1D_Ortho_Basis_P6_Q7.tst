<?xml version="1.0" encoding="utf-8"?>
<test>
    <description>Project1D and Diff Segment Orthogonal basis P=6 Q=7</description>
    <executable>LocProject</executable>
    <parameters>-s Segment -b Ortho_A -o 6 -p 7 -c 0.0 1.0 -d</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-12">9.19829e-15</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-12">3.90799e-14</value>
        </metric>
    </metrics>
</test>


