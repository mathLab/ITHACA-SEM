<?xml version="1.0" encoding="utf-8" ?>
<test>
    <description>StdProject3D Pyramid Modified basis P=6 Q=7</description>
    <executable>StdProject</executable>
    <parameters>-s pyramid -b Modified_A Modified_A ModifiedPyr_C -o 6 6 6 -p 7 7 6</parameters>
    <metrics>
        <metric type="L2" id="1">
            <value tolerance="1e-9">3.2756e-12</value>
        </metric>
        <metric type="Linf" id="2">
            <value tolerance="1e-9">9.72786e-13</value>
        </metric>
    </metrics>
</test>
